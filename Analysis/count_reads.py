
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
from helper_functions import unique, make_regexes, make_sgIDDict, getsgID, hamming_distance, revcom, read_project_info, read_parameter_info

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

project_info_file = root + "/Parameters/" + project_name + "_project_file.txt"
parameter_info_file = root + "/Parameters/" + parameter_name + "_parameter_file.txt"
project_info = read_project_info(project_info_file)
parameter_info = read_parameter_info(parameter_info_file)

print("In count_reads, sample is {}\n".format(sample))
sgIDDict = make_sgIDDict(project_info.sgids, project_info.sgRNAs) #keys are sequences, values are genes targeted
sgIDDict_invert = {value : key for (key, value) in sgIDDict.items()}

########################################################################################
# 3. Set up output Files
########################################################################################


naming_stub = project_name + "_" + parameter_name + "_" + sample
record_name = root + "/" + project_name + "/" + parameter_name + "/records/" + naming_stub + "_record.txt"
reject_read_name = root + "/" + project_name + "/" + parameter_name + "/rejected/" + naming_stub + "_rejects.txt"
output_name = root + "/" + project_name + "/" + parameter_name + "/raw_counts/" + naming_stub + "_counts.txt"
record = open(record_name, 'wt')
reject = open(reject_read_name, 'wt')
output = open(output_name, 'wt')



reject.write("Reason,line_seq1,line_seq2\n")

########################################################################################
# 4. Set up regular expressions for testing
########################################################################################

regexes_to_check = make_regexes(parameter_info.R1_regex_looseness,parameter_info.R2_regex_looseness)
print("In count_reads, the regular expressions I'm testing are {}\n".format(regexes_to_check))


#########################################################################################
# 5. main loop, running through lines of file
########################################################################################

File1 = sample_path
print("File1 is {}\n".format(File1))

if "_R1" in File1:
    File2 =File1.replace('_R1','_R2')
elif "_1.fq" in File1:
    File2 = File1.replace('_1.fq','_2.fq')
else:
    sys.exit("Problem, don't have way to differentiate forward and reverse read files\n")

record.write("Tuba seq pipeline\n")
record.write("Beginning Part 1 of Tuba seq pipeline: counting reads\n ")
record.write("Timestamp: {}\n".format(now.strftime("%Y-%m-%d %H:%M")))
record.write("Project: {}, Parameters: {}\n".format(project_name, parameter_name))
record.write("Sample: {},\n\n".format(sample))
record.write("Input files:\nRead1: {}\nRead2: {}\n\n".format(File1, File2))
record.write("Output files:\nPrimary: {}\nRejected Reads: {}\n Record: {}\n".format(output_name,reject_read_name, record_name))
record.write("Regexes being searched for in each read:\n")
record.write("R1: {}, R2: {}\n".format(regexes_to_check[0], regexes_to_check[1]))

f1 = gzip.open(File1,'rt')
f2 = gzip.open(File2,'rt')

None_seq_dict = {}
sgIDBCdict = {}
line_seq1 = 1

total_reads = 0
n_tumors = 0
accepted_reads = 0
rejected_reads = 0

R1_didnt_match = 0
R2_didnt_match = 0
regex_rejects = 0
BC_match = 0
BC_match_fail = 0

prelim_none_counts = 0
prelim_good_counts = 0


while (line_seq1):
    line1 = f1.readline().rstrip() # skip the first line
    print("line 1 is {}\n".format(line1))

    line_seq1 = f1.readline().rstrip() 
    line1 = f1.readline().rstrip() # skip the third line
    line_qua1 = f1.readline().rstrip() # get the sequencing quality
    line2 = f2.readline().rstrip()
    line_seq2 = f2.readline().rstrip()
    line2 = f2.readline().rstrip()
    line_qua2 = f2.readline().rstrip()
  #  if project_info.company in ["Novogene","NovoGene","novogene","NOVOGENE","NG"]:
 #       print("checking that I am responding to company cue\n")
 #       line_seq2_comp = line_seq2
 #   else:
 #       line_seq2_comp=revcom(line_seq2)
    total_reads+=1
    match = 0

    regexR1 = re.compile(regexes_to_check[0])
    regexR2 = re.compile(regexes_to_check[1])
    #    if regexR1.search(line_seq1) and regexR2.search(line_seq2_comp):

    if regexR1.search(line_seq1) and regexR2.search(line_seq2):
        match += 1
        k1 = regexR1.search(line_seq1) # align R1
        k2 = regexR2.search(line_seq2) # align R2
    elif regexR1.search(line_seq1) and regexR2.search(revcom(line_seq2)):
        match += 1
        k1 = regexR1.search(line_seq1) # align R1
        k2 = regexR2.search(revcom(line_seq2)) # align R2  
    else:
        rejected_reads += 1
        regex_rejects += 1
        reject.write("R1orR2_didnt_match,{},{}\n".format(line_seq1,line_seq2))

    if match ==1:
        R1sgID = k1.group(1) # upstream of R1
        R1BC = k1.group(2) # downstream of R1
        sgID_1 = getsgID(sgIDDict, R1sgID, 2) #note to self: getsgID will assign none to anything that is not in the expected sgids for this project.
        R2sgID = k2.group(1) # upstream of R1
        R2BC = k2.group(2) # downstream of R2_RC
        sgID_2 = getsgID(sgIDDict, R2sgID, 2)
        
        BC_dist=hamming_distance(R1BC,R2BC)
        #print("R1BC is {}, R2BC is {}, R1sgid is {}, r2sgid is {}, BC_dist is {}\n".format(R1BC, R2BC, sgID_1, sgID_2, BC_dist))

        if BC_dist <= int(parameter_info.dist_between_bc_for_read) and ('N' not in R1BC): #barcodes match
            BC_match += 1
            myKey = sgID_1 + "," + R1BC
            accepted_reads += 1
            if myKey in sgIDBCdict:
                sgIDBCdict[myKey] += 1
            else:
                n_tumors += 1
                sgIDBCdict[myKey] = 1

            if sgID_1 == "None":
                None_seq_dict[myKey]=R1sgID

        else:
            BC_match_fail += 1
            rejected_reads += 1
            reject.write("BC_match_fail,{},{}\n".format(line_seq1,line_seq2))


f1.close()
end_time = time.monotonic()
td = timedelta(seconds=end_time - start_time)
percent_accepted_including_none = float(accepted_reads)/float(total_reads)
percent_rejected = float(rejected_reads)/float(total_reads)


record.write("Counting reads took: {}\n\n".format(td))
record.write("##########\n")
record.write("Reads processed: {}. Tumors processed {}\n".format(total_reads, n_tumors))
record.write("Reads accepted including those in tumors with unexpected sgID sequence: {} ({:.2%})\n".format(accepted_reads,percent_accepted_including_none ))
record.write("Reads rejected: {} ({:.2%})\n".format(rejected_reads, percent_rejected))

record.write("Reads discarded at regex matching stage: {}\n".format(regex_rejects))
record.write("Reads discarded at barcode matching stage: {}\n".format(BC_match_fail))
record.write("Reads where R1-R2 barcodes match: {}\n".format(BC_match))
record.write("##########\n")

output.write("sgID,BC,Count,sgIDseq\n")
for k,v in sorted(sgIDBCdict.items()):
    if k.strip().split(",")[0]=="None":
        prelim_none_counts+=int(v)
        output.write("{},{},{}\n".format(k,v,None_seq_dict[k]))
    else:
        prelim_good_counts+=int(v)
        output.write("{},{},{}\n".format(k,v,sgIDDict_invert[k.strip().split(",")[0]]))
record.write("Prelim None read count: {}\n".format(prelim_none_counts))
record.write("Prelim good read count: {}\n".format(prelim_good_counts))
