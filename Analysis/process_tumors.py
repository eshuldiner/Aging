'''
This script reads in all the tumors from the final processed files into one table and annotates 
with useful information. 

It creates jitterplots while doing this

'''

#python3 process_tumors.py --project_name=UCSF_Injury_corr2 --parameter_name=2 --root=/labs/mwinslow/Emily/

from helper_functions import read_project_info, read_parameter_info, get_gc, pull_sgid, pull_bc
from barcode_diversity_functions import barcodeDiversityCutoffs, identifyKeepers, removeContamination
import getopt
import glob
import sys
import statistics
from subprocess import call
import numpy as np

####################################################
# Step 1, Take arguments, get project info
####################################################

try:
    opts, args = getopt.getopt(sys.argv[1:], "", ["root=", 
        "project_name=", 
        "parameter_name="])
except getopt.GetoptError:
    print("no arguments recognized\n")
    sys.exit(2)

for opt, arg in opts:
    if opt in ("--root"):
        root = arg
        print("found root, {}\n".format(root))
    elif opt in ("--project_name"):
        project_id = arg
        print("found project_name, {}\n".format(project_id))
    elif opt in ("--parameter_name"):
        parameter_name = arg
        print("found parameter_name, {}\n".format(parameter_name))
    else:
        assert False, "unhandled option"

project_info_file = root + "/Parameters/" + project_id + "_project_file.txt"
parameter_info_file = root + "/Parameters/" + parameter_name + "_parameter_file.txt"
outfname_no_cutoff = root + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_tumors_unfiltered.txt"
outfname_no_cutoff_filtered = root + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_tumors.txt"
contam_report_f= root + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_contamination_report.txt"
spi_out_name = root + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_spike_ins.txt"

o_all = open(outfname_no_cutoff, 'wt')
s_o = open(spi_out_name, 'wt')

r_jitter_script = root + "/Visualization/jitter_plotter.R"
contamination_plotter = root + "/Visualization/contamination_plotter.R"


project_info = read_project_info(project_info_file)
parameter_info = read_parameter_info(parameter_info_file)
tumor_output_dir = root + project_id + "/" + parameter_name + "/filtered/"

####Creating a single tumor file that contains all the tumors across all samples
tumor_data = {}
for_bc_diversity_analysis = {}
nsamples=0

o_all.write("sgID,BC,Count,CellNum,Gene,Sample,Genotype,Titer,GC,Treatment\n")
s_o.write("sgID,BC,Count,CellNum,Gene,Sample,Genotype,Titer,GC,Treatment\n")
for s in project_info.sample_ids:

	nsamples +=1
	print("Processing data from sample {}, {}\n".format(s, project_info.sample_to_treatment[s]))
	fname = tumor_output_dir + project_id + "_" + parameter_name + "_" + s + "_final.txt"
	f = open(fname,'rt')
	f.readline()
	for l in f:
		fields = l.strip().split(",")
		if len(fields) >1:
			sgid = fields[0]
			bc = fields[1]
			rc = int(fields[2])
			cellnum = float(fields[3])
			gc = get_gc(bc)
			sgid_bc = sgid + "_" + bc

			if sgid in for_bc_diversity_analysis:
				for_bc_diversity_analysis[sgid].append(bc)
			else:
				for_bc_diversity_analysis[sgid] = [bc]

			if sgid_bc in tumor_data:
				tumor_data[sgid_bc]["SAMPLES"].append(s)
				tumor_data[sgid_bc]["RC"].append(rc)
			else:
				tumor_data[sgid_bc] = {}
				tumor_data[sgid_bc]["SAMPLES"] = [s]
				tumor_data[sgid_bc]["RC"] = [rc]
			o_all.write("{},{},{},{},{},{},{}\n".format(l.strip(), project_info.sgids_to_gene[sgid], s, project_info.sample_to_gt[s], project_info.sample_to_titer[s], gc, project_info.sample_to_treatment[s]))
			if sgid == "Spi":
				s_o.write("{},{},{},{},{},{},{}\n".format(l.strip(), project_info.sgids_to_gene[sgid], s, project_info.sample_to_gt[s], project_info.sample_to_titer[s], gc,project_info.sample_to_treatment[s]))
o_all.close()


#############################################
##########Contamination correction###########
#############################################

cutoffs = barcodeDiversityCutoffs(root, project_id, parameter_name, float(parameter_info.contam_removal_threshold), for_bc_diversity_analysis, nsamples)
keepers = identifyKeepers(tumor_data)
removeContamination(outfname_no_cutoff, outfname_no_cutoff_filtered, cutoffs, contamination_report=contam_report_f, keepers=keepers)

#Plots
call("module load r; R --vanilla --args " + contam_report_f + " " + str(project_id) + " " + str(parameter_name) + " " + root +  " < " + contamination_plotter, shell = True)
call("module load r; R --vanilla --args " + outfname_no_cutoff_filtered +  " " + str(",".join(project_info.inerts)) + " " + root + " " + project_id + " " + parameter_name + " < " + r_jitter_script, shell = True)
