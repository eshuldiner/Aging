'''
This script harvest information from the records that are written during read processing
'''

#python3 analyze_records.py --project_name=UCSF_Injury_corr3 --parameter_name=2 --root=/labs/mwinslow/Emily/

from helper_functions import read_project_info, read_parameter_info
import getopt
import glob
import sys
from subprocess import call

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

project_info_file = root + "/tubaseq_project_files/" + project_id + "_project_file.txt"
parameter_info_file = root + "/tubaseq_parameter_files/" + parameter_name + "_parameter_file.txt"
project_info = read_project_info(project_info_file)
parameter_info = read_parameter_info(parameter_info_file)
r_record_script = root + "/Visualization/record_plotter.R"

record_dir = root + "/" +project_id + "/" + parameter_name + "/records/"

# Dictionary of fields I am actually interested in; values are strings to search in each line

sample_info_dict = {}
sample_index_storage = {}

of_interest = {
			"Processed" : "Reads processed", 
			"Accepted_with_None" : "Reads accepted including those in tumors with unexpected sgID sequence",
			"Rejected_read_stage" : "Reads rejected", 
			"Reject_regex" : "discarded at regex matching stage", #good
			"Reject_BCmatch" : "discarded at barcode matching stage",
			"Pass_BCmatch" : "Reads where R1-R2 barcodes match",	
			"Other_project_reads":"Reads mapping to other projects",
			"Insert_Length_Reject":"Reads discarded due to incorrect insert length",
			"Not_other_project_reads":"Reads NOT mapping to other projects",
			"None_reads" : "Reads not mapping to this project",
			"Accepted_no_None" : "Reads excluding none/other projects",
			}


for s in project_info.sample_ids:
	try:
		fname = record_dir + project_id + "_" + parameter_name + "_" + s + "_record.txt"

		values_of_interest = {
				"Processed" : 0, #good
				"Accepted_with_None" : 0,
				"Rejected_read_stage" : 0, #good
				"Reject_regex" : 0, #good
				"Reject_BCmatch" : 0,
				"Pass_BCmatch" : 0,	
				"Other_project_reads":0,
				"Insert_Length_Reject":0,
				"Not_other_project_reads":0,
				"None_reads":0,
				"Accepted_no_None": 0,
					}
		f = open(fname,'rt')

		for l in f:
			print("l is {}\n".format(l))
			fields = l.strip().split(":")
			print("fields is {}\n".format(l))
			for val, searcher in of_interest.items():
				print("working on val, searcher {}\n".format(val, searcher))
				if searcher in fields[0]:
					print('finding searcher')
					values_of_interest[val] = fields[1].strip().strip(".")


		sample_info_dict[s] = values_of_interest
	except:
		print("problem with sample {}\n".format(s))
		pass
#Now you have dictionary containing info collected in the record file for each sample. You should write this to a table

include = [	"Processed", "Accepted_with_None","Accepted_no_None", "Rejected_read_stage", "Reject_regex", "Insert_Length_Reject","Reject_BCmatch", "Pass_BCmatch","None_reads","Other_project_reads","Not_other_project_reads"]
record_summary_out_name = root + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_record_summary.txt"
record_summary_out = open(record_summary_out_name, 'wt')
record_summary_out.write("Sample,GT")
for x in include:
	record_summary_out.write(",{}".format(x))
record_summary_out.write("\n")

first = True
for sample in project_info.sample_ids:
	if first == False:
		record_summary_out.write("\n")
	first=False
	record_summary_out.write("{},{}".format(sample,project_info.sample_to_gt[sample]))
	for x in include:
		record_summary_out.write(",{}".format(sample_info_dict[sample][x]))
	
record_summary_out.close()
call("module load r; R --vanilla --args " + record_summary_out_name + " " + str(project_id) + " " + str(parameter_name) + " " + str(",".join(project_info.inerts)) + " " + root +  " < " + r_record_script, shell = True)





