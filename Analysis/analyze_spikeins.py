'''
This script harvest information from the records that are written during read processing
'''

#python3 analyze_spikeins.py --project_name=UCSF_Injury --parameter_name=2 --root=/labs/mwinslow/Emily/

from helper_functions import read_project_info, read_parameter_info
import getopt
import glob
import sys
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

project_info_file = root + "/tubaseq_project_files/" + project_id + "_project_file.txt"
parameter_info_file = root + "/tubaseq_parameter_files/" + parameter_name + "_parameter_file.txt"
project_info = read_project_info(project_info_file)
parameter_info = read_parameter_info(parameter_info_file)
sample_to_depth = project_info.sample_to_depth
r_spi_plot_script = root + "/tubaseq_pipeline/R_post_processing/spi_plotter.R"

print("In analyze spike-ins")

infname = root + "/" + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_spike_ins.txt"
count_outfname = root + "/" + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_spike_in_count.txt"
top_bc_outfname = root + "/" + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_spike_in_top_bcs.txt"
data = open(infname,'rt')
header=data.readline().strip().split(",")

n_out = open(count_outfname, 'wt')
n_out.write("Sample,Genotype,Treatment,totalSpikeInReads,NSpikeIns,N_spike_gt_10,N_spike_gt_50,depth,intended_rel_depth\n")

top_bc_out = open(top_bc_outfname, 'wt')
top_bc_out.write("Sample,Genotype,Barcode,Reads,percentReads\n")


total_bc_reads = {}
n_spikes_per_sample = {}
n_spikes_per_sample_gt10 = {}
n_spikes_per_sample_gt50 = {}

unique_spi_bc = {}
unique_spi_bc_gt10 = {}
unique_spi_bc_gt50 = {}
depth = {}

unique_spi_rc = {}
unique_spi_rc_gt10 = {}
unique_spi_rc_gt50 = {}


print("LOOK FOR ME!\n")
for l in data:
	fields = l.strip().split(",")
	print("line is: {}\n".format(l))
	bc = fields[header.index("BC")]
	rc = int(fields[header.index("Count")])
	cn =float(fields[header.index("CellNum")])
	sample = fields[header.index("Sample")]
	depth[sample] = rc/cn
	if sample in n_spikes_per_sample.keys(): #you've seen this sample before
		n_spikes_per_sample[sample] += 1
		unique_spi_bc[sample].append(bc)
		unique_spi_rc[sample].append(rc)
		total_bc_reads[sample]+=rc
	else:
		n_spikes_per_sample[sample] = 1
		unique_spi_bc[sample] = [bc]
		unique_spi_rc[sample] = [rc]
		total_bc_reads[sample] = rc
	if rc >10:
		if sample in n_spikes_per_sample_gt10.keys(): #you've seen this sample before
			n_spikes_per_sample_gt10[sample] += 1
			unique_spi_bc_gt10[sample].append(bc)
			unique_spi_rc_gt10[sample].append(rc)
		else:
			n_spikes_per_sample_gt10[sample] = 1
			unique_spi_bc_gt10[sample] = [bc]
			unique_spi_rc_gt10[sample]= [rc]
	if rc >50:
		if sample in n_spikes_per_sample_gt50.keys(): #you've seen this sample before
			n_spikes_per_sample_gt50[sample] += 1
			unique_spi_bc_gt50[sample].append(bc)
			unique_spi_rc_gt50[sample].append(rc)
		else:
			n_spikes_per_sample_gt50[sample] = 1
			unique_spi_bc_gt50[sample] = [bc]
			unique_spi_rc_gt50[sample] = [rc]


for sample in project_info.sample_ids:
	if sample not in n_spikes_per_sample.keys():
		n_spikes_per_sample[sample] = 0
	if sample not in n_spikes_per_sample_gt10.keys():
		n_spikes_per_sample_gt10[sample] = 0
	if sample not in n_spikes_per_sample_gt50.keys():
		n_spikes_per_sample_gt50[sample] = 0

for sample in n_spikes_per_sample.keys():
	try:
		n_out.write("{},{},{},{},{},{},{},{},{}\n".format(sample, project_info.sample_to_gt[sample], project_info.sample_to_treatment[sample], total_bc_reads[sample],n_spikes_per_sample[sample],n_spikes_per_sample_gt10[sample],n_spikes_per_sample_gt50[sample], depth[sample], sample_to_depth[sample]))
	except:
		print("Something went wrong for sample {}\n".format(sample))

def determine_bc_to_plot(intended, n_spikes_per_sample):
	n_spi_all = []
	for sample in n_spikes_per_sample.keys():
		n_spi_all.append(n_spikes_per_sample[sample])
	print("n_spi_all is {}\n".format(n_spi_all))
	print("max is: {}\n".format(max(n_spi_all)))
	bc_to_plot = min(intended, max(n_spi_all))
	return(bc_to_plot)

n_bc = determine_bc_to_plot(10, n_spikes_per_sample)

#unique_spi bc and rc are in descending rc order

tracking_bc = []
trimmed_top_bc = {}
trimmed_top_rc = {}

print("n_bc is {}\n".format(n_bc))

for sample in unique_spi_rc_gt50.keys():
	top_bc = unique_spi_bc_gt50[sample][0:n_bc]
	top_rc = unique_spi_rc_gt50[sample][0:n_bc]
	tracking_bc.extend(top_bc)
	trimmed_top_bc[sample] = top_bc
	trimmed_top_rc[sample] = top_rc

print("BEFORE FILLING IN WITH ZEROS!")
for sample in trimmed_top_bc.keys():
	print("sample is {}, bc {}, rc {}\n".format(sample, trimmed_top_bc[sample], trimmed_top_rc[sample]))


###Now make sure that I have a measurement for every sample for all the top bc
for bc in set(tracking_bc):
	for sample in trimmed_top_rc.keys():
		if bc not in trimmed_top_bc[sample]:
			trimmed_top_bc[sample].append(bc)
			trimmed_top_rc[sample].append(0)

print("AFTER FILLING IN WITH ZEROS!")
for sample in trimmed_top_bc.keys():
	print("sample is {}, bc {}, rc {}\n".format(sample, trimmed_top_bc[sample], trimmed_top_rc[sample]))


reformat_info = {}

for sample in trimmed_top_bc.keys():
	reformat_info[sample] = {}
	for i, bc in enumerate(trimmed_top_bc[sample]):
		reformat_info[sample][bc] = trimmed_top_rc[sample][i]


for sample in reformat_info.keys():
	for bc in reformat_info[sample].keys():
		#print("{},{},{}\n".format(sample,bc, reformat_info[sample][bc]))
		top_bc_out.write("{},{},{},{},{}\n".format(sample, project_info.sample_to_gt[sample],bc, reformat_info[sample][bc], reformat_info[sample][bc]/total_bc_reads[sample]))




n_out.close()
top_bc_out.close()

genotypes_in_order = []

for s in project_info.sample_ids:
    genotypes_in_order.append(project_info.sample_to_gt[s])
call("module load r; R --vanilla --args " + top_bc_outfname + " " + count_outfname + " " + str(project_id) + " " + str(parameter_name) + " " + str(",".join(genotypes_in_order)) + " " + root +  " < " + r_spi_plot_script, shell = True)




