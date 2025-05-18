########################################################################################
# 1. Load packages
########################################################################################
import sys
import os
import regex as re
import gzip
import getopt
import time
import datetime
import numpy as np
from datetime import timedelta
from subprocess import call
from helper_functions import hamming_distance, read_project_info, ListToDict, onlyGoodSpikeIns, get_conversion_factor_unweighted
import numpy as np


now = datetime.datetime.now()
start_time= time.monotonic()

########################################################################################
# 2. Take arguments and build necessary resources
########################################################################################
try:
    opts, args = getopt.getopt(sys.argv[1:], "s:", ["sample=", 
        "sample_path=", 
        "root=", 
        "project_name=", 
        "parameter_name="])
except getopt.GetoptError:
    print("no arguments recognized\n")
    sys.exit(2)

for opt, arg in opts:
    if opt in ("--sample"):
        sample = arg
    elif opt in ("--sample_path"):
        sample_path = arg
    elif opt in ("--root"):
        root = arg
    elif opt in ("--project_name"):
        project_name = arg
    elif opt in ("--parameter_name"):
        parameter_name = arg
    else:
        assert False, "unhandled option"


#get all the info for this project:
project_info_file = root + "/Parameters/" + project_name + "_project_file.txt"
project_info = read_project_info(project_info_file) 
########################################################################################
# 3. Set up input/output Files
########################################################################################

naming_stub = project_name + "_" + parameter_name + "_" + sample
input_name = root + "/" + project_name + "/" + parameter_name + "/filtered/" + naming_stub + "_clustered.txt"
record_name = root + "/" + project_name + "/" + parameter_name + "/records/" + naming_stub + "_record.txt"
output_name = root + project_name + "/" + parameter_name + "/filtered/" + naming_stub + "_final.txt"
record = open(record_name, 'a+')
output = open(output_name, 'wt')
_input = open(input_name, 'rt')
_input.readline()

output.write("sgid,bc,rc,cellnum\n")

caution = False
errors = []
record.write("In convert_to_cells.py\n")


#########################################################################################
# 4. read in data
########################################################################################

'''
Note to self: tumors_all here is keyed sgid --> sgid,bc,rc (as a string)
'''
tumors_in = 0
tumors_all = {}

for l in _input:
    fields = l.strip().split(",")
    sgid = fields[0]
    bc = fields[1]
    rc = fields[2]
    tumors_in += 1
    key = sgid + "," + bc + "," + rc
    if sgid in tumors_all.keys():
        tumors_all[sgid].append(key)
    else:
        tumors_all[sgid] = [key]

record.write("Read in data for Sample {}\n".format(sample))
#########################################################################################
# 5. Analyze the spike-ins (get read count --> cell number conversion factor)
########################################################################################
expected_spi_bc = project_info.sample_to_spike_bcs[sample]
spi_bc_sizes = project_info.sample_to_spike_sizes[sample]

# Checks that the spike-in information provided by the user is internally consistent.
if len(expected_spi_bc) != len(spi_bc_sizes):
    record.write("ERROR WITH SPIKE-INS! you expect {} spike-ins, but only have sizes for {}\n".format(len(expected_spi_bc), len(spi_bc_sizes)))
    print("ERROR WITH SPIKE-INS! you expect {} spike-ins, but only have sizes for {}\n".format(len(expected_spi_bc), len(spi_bc_sizes)))
    sys.exit(2)

#Spike-in barcodes will automatically be in descending read count order

#Check that the number of spike-ins associated with the sample is at least as large as the expected # of spike-in barcodes:
bc_to_sizes_theoretical = dict(zip(expected_spi_bc, spi_bc_sizes))
if "Spi" in tumors_all:
    record.write("Successfully identified spike-ins\n")
    SpikeInData = tumors_all["Spi"]
    bc_to_size = ListToDict(SpikeInData)
    if (len(SpikeInData) < len(spi_bc_sizes)): #there are spike-ins, but fewer than the expected #
        print("The total number of spike-in barcodes ({})is less than the expected number ({}) \n".format(str(len(SpikeInData)),str(len(spi_bc_sizes)))) 
        caution = True
        errors.append("Total # of barcodes with spike-in sgID < expected number")
        if onlyGoodSpikeIns(bc_to_size, expected_spi_bc): #Not enough bc but they are all bc you expect
            #do calculation using RC from all spike-ins in SpikeInData
            factor = get_conversion_factor_unweighted(bc_to_sizes_theoretical, bc_to_size)
        else: #You don't have enough bc, and some of them are unexpected bc
            print("Some of the spike-in barcodes are not those that you expected\n")
            #do calculation using RC from the spike-ins that you have
            to_use = { bc:rc for (bc,rc) in bc_to_size.items() if bc in expected_spi_bc }
            factor = get_conversion_factor_unweighted(bc_to_sizes_theoretical, to_use)
            errors.append("Unexpected spike-in barcodes are more highly represented than real spike-in barcodes")
    else: #there are at least the expected number of barcodes
        print("The total number of spike-in barcodes is at least as large as the expected number\n")
        largest_spike_ins = []
        for i in range(len(spi_bc_sizes)):
            largest_spike_ins.append(SpikeInData[i])
        if onlyGoodSpikeIns(largest_spike_ins, expected_spi_bc):
            print("The X most highly represented spike-in barcodes are as expected\n")
            caution = False
            #complete success case: the 3 largest spike-ins are the expected barcodes.
            to_use = { bc:rc for (bc,rc) in bc_to_size.items() if bc in expected_spi_bc } #used to do order-based filtering, now doing it explicitly

            ##Section added April 26 2023 to flag samples with extreme outlier spike-ins:
            quick_cells_per_read = {}
            rel_to_median= {}
            for bc in to_use:
                quick_cells_per_read[bc]=bc_to_sizes_theoretical[bc]/int(to_use[bc])
            median = np.median(list(quick_cells_per_read.values()))

            for bc in to_use:
                rel_to_median[bc]=quick_cells_per_read[bc]/median

            to_use_filtered = {}
            problem_counter = 0
            for bc,val in rel_to_median.items():
                if val > 10 or val < 0.1: # if a spike-in barcode is more than 10x away from the median, take it out.
                    problem_counter+=1
                    cation=True
                else:
                    to_use_filtered[bc] = to_use[bc]

            factor = get_conversion_factor_unweighted(bc_to_sizes_theoretical, to_use_filtered)

        else:
            print("Some of the most highly represented spike-in barcodes are not those that you would expect\n")
            caution = True
            if all(item in bc_to_size  for item  in expected_spi_bc): #They are not the largest 3, but you do see all expected barcodes
                print("All X expected barcodes are present in the dataset\n")
                errors.append("ALL SPIKE-IN BC PRESENT BUT AT LEAST ONE UNEXPECTED BC IS MORE HIGHLY REPRESENTED")
                to_use = { bc:rc for (bc,rc) in bc_to_size.items() if bc in expected_spi_bc }
                factor = get_conversion_factor_unweighted(bc_to_sizes_theoretical, to_use)

                #do calculation using all expected spike-in barcodes
            else: #You do not see all N expected barcodes
                print("Not all expected barcodes are present\n")
                errors.append("NOT ALL SPIKE-IN BARCODES PRESENT")
                to_use = { bc:rc for (bc,rc) in bc_to_size.items() if bc in expected_spi_bc }
                factor = get_conversion_factor_unweighted(bc_to_sizes_theoretical, to_use)
                #do calculation RC from the spike-ins that you have 

else:
    #there were no spike-in bc associated with this sample
    print("No spike-in barcodes present\n")
    caution = True
    errors.append("No spike-ins associated with sample; using read count of 1 for conversion")
    factor = sum(spi_bc_sizes)/1


#########################################################################################
# 6. Re-write file with cell number
########################################################################################
for sgid in tumors_all.keys():
    for tumor in tumors_all[sgid]:
        fields=tumor.strip().split(",")
        cellnum = int(fields[2])*factor
        output.write("{},{}\n".format(tumor,cellnum))

#########################################################################################
# 7. Write out-come of RC -> cell number conversion to record file
#########################################################################################
#print(f"{'Name' : <10}{'Marks' : ^10}{'Division' : ^10}{'ID' : >5}")
record.write("SAMPLE,STATUS,ERRORS\n")
if caution:
    record.write("{},CAUTION,{}\n".format(sample, ';'.join(errors)))
else:
    record.write("{},SUCCESS,None\n".format(sample))

