'''
This script estimates a lambda and an overall barcode diversity (D) for each sgID in a Tuba-seq dataset by assuming that the probability of each BC occuring in i mice is poisson-distributed

How this works:

	Can calculate directly from data:
		-mu_nonzero = the average # of occurences of each BC (by necessity, across mice)
		-L = # of draws from the pool of possible barcodes, we are estimating as the total # of BC associated with the sgID

	Then use mu_nonzero to estimate lambda for the poisson, and use lambda to estimate D based on lambda=L/D

Given lambda, can then calculate the prob of seing a BC in >X samples, and use that to decide on a cut-off for when to 
discard BC that are present in too many samples.

You can also use lambda to calculate the probability of seeing the same barcode more than once in a sample.

'''
import numpy as np
import sys
from helper_functions import pull_sgid, pull_bc


def input_data_for_bc_diversity_analysis(tumor_file_name):
	input_data=open(tumor_file_name,'rt')
	input_data.readline()

	data = {}
	all_samples = []

	for l in input_data:
		fields = l.strip().split(",")
		sgID = fields[0]
		BC = fields[1]
		sample=fields[5]
		if sample not in all_samples:
			all_samples.append(sample)

		if sgID in data:
			data[sgID].append(BC)
		else:
			data[sgID]=[BC]
	input_data.close()
	return(data,len(all_samples))




def get_mu_nonzero(data, sgID):
	'''calculate the average number of occurences of each BC (across mice) for an sgID.
	This is done by calculating the # of mice each BC occur in, and then taking the average'''
	bc_counts = {}

	bc = data[sgID]

	for b in bc:
		if b in bc_counts:
			bc_counts[b]+=1
		else:
			bc_counts[b]=1

	counts = list(bc_counts.values())
	mu_nonzero = np.mean(counts)

	return(mu_nonzero)



def get_L(data, sgID):
	''' Returns the # of tumors associated with an sgID (multi-counting duplicate BC)'''
	return(len(data[sgID]))

def get_lambda(mu_nonzero,N_samp):
	'''Given mu_nonzero, calculate lambda from equation mu_nonzero = lambda/(1-e^-lambda).
	This will need to be down empirically (i.e. with an optimizer)'''
	from scipy.optimize import brentq
	def fun(_lambda):
		return(_lambda-mu_nonzero+mu_nonzero*np.exp(-1*_lambda))

	lambda_sol = brentq(fun, 0.0000001, N_samp) # need to start search above 0 because 0 is also a root
	return(lambda_sol)

def get_D(L,_lambda):
	return(L/_lambda)

def get_cutoff(certainty, _lambda, nsamples):
	'''Given the poisson distribution defined by _lambda, and the specified certainty, how many
	samples is too many samples for a BC to occur in?

	e.g. if certainty = 1% then what is the number of samples such that given the poisson,
	99% of BC should be in fewer samples?

	Reminder that for the poisson cdf the probability returned  includes the limit, i.e. is to see the limit OR less
	'''
	from scipy.stats import poisson
	for i in range(nsamples+1):
		cumulative_prob = poisson.cdf(k=i,mu=_lambda)
		if 1-cumulative_prob<certainty:
			return(i+1)


def barcodeDiversityCutoffs(root, project_id, parameter_name, certainty, for_bc_diversity_analysis, nsamples):
	##PER SGID INFO, FOCUSED ON BARCODE DIVERSITY
	bc_diversity_outf = root + project_id + "/" + parameter_name + "/summary/" + project_id + "_" + parameter_name + "_Barcode_Diversity_Report.txt"
	bc_diversity_output = open(bc_diversity_outf, 'wt')
	bc_diversity_output.write("sgID,mu_nz,L,D,lambda,N_samp_in_exp,CutoffCalc,Cutoff\n")

	cutoffs = {}

	for sgid in for_bc_diversity_analysis.keys():
		if sgid !="Spi":
			mu_nz = get_mu_nonzero(for_bc_diversity_analysis, sgid) #average number of occurences of each BC (across mice) for an sgID.
			L = get_L(for_bc_diversity_analysis, sgid) #the # of tumors associated with an sgID (multi-counting duplicate BC)
			try:
				_lambda = get_lambda(mu_nz, nsamples) #lambda is I think the mean of the poisson
				D=get_D(L, _lambda) #I think D is supposed to be the estimated number of unique barcodes. L/lambda
				cutoff_calc = get_cutoff(certainty ,_lambda, nsamples)
				cutoff = min(cutoff_calc,99)

				cutoffs[sgid] = cutoff
			except:
				_lambda = "NA"
				D = "NA"
				cutoff_calc = 99
				cutoff = 99
				cutoffs[sgid] = 99
			bc_diversity_output.write("{},{},{},{},{},{},{},{}\n".format(sgid, mu_nz, L, D, _lambda,nsamples, cutoff_calc,cutoff ))

	cutoffs["Spi"]=10000 #setting cutoff arbitrarily high for Spike-ins so that they are never removed based on being in too many samples

	bc_diversity_output.close()
	return(cutoffs)

def identifyKeepers(tumor_data):
	counter=0
	perLargest_keepers = {}

	for s in tumor_data:
	
		splitter = s.strip().split("_")
		samples = np.array(tumor_data[s]["SAMPLES"])
		sizes = np.array(tumor_data[s]["RC"])
		#sgid = splitter[0] + "_" + splitter[1]
		sgid = "_".join(splitter[:-1])

		if len(sizes)>1:
			counter+=1

			maxSize = sizes.max()
			totalReads = np.sum(sizes)
			numSamples = len(samples)
			variance = np.var(sizes)

			sizes_adj = sizes/totalReads
			perLargest=sizes_adj.max()

			if perLargest > 0.95:
				perLargest_keepers[s] = list(samples)[list(sizes_adj).index(perLargest)]

			if counter%1000==0:
				print("{} barcodes processed\n".format(counter))
			return(perLargest_keepers)


def removeContamination(infilename, outfilename, removal_threshold_dict, contamination_report, keepers):
	'''
	Remove barcodes that occur in "too many" samples and are therefore likely to be contamination
	(based on analysis of barcode diversity). But if  >95% of  reads are assigned to a single sample, 
	then retain the barcode in the sample with the bulk of the sequencing reads and discard in all other samples. 
	'''

	infile = open(infilename,'rt')
	counter = 0
	tumors_in = {}
	header=infile.readline()
	#tumors_out.write(header)
	for l in infile:
		counter+=1
		fields=l.strip().split(",")
		sgid_bc = fields[0] + "_" + fields[1]
		if sgid_bc in tumors_in:
			tumors_in[sgid_bc].append(l.strip())
		else:
			tumors_in[sgid_bc]=[l.strip()]
	infile.close()

	# step 2: write or don't tumors based on # of occurences. Keep track of it.
	outfile=open(outfilename,'wt')
	if contamination_report !="NA":
		contam_report_out=open(contamination_report,'wt')
		contam_report_out.write("sgID_BC,sgID,Cutoff_for_sgID,BC,N_samples,N_reads,Rejected\n")
	
	outfile.write("{}".format(header))
	for sgid_bc in tumors_in: #egs 8/18/2021: modifying code so that all sgid-bc that occur in more than one sample are written and it's annotated whether they are removed or not based on the cutoff
		reads = 0
		sgid = pull_sgid(sgid_bc)
		removal_threshold = removal_threshold_dict[sgid]
		for tl in tumors_in[sgid_bc]:
			reads+=int(tl.split(",")[2])
		if len(tumors_in[sgid_bc])>=removal_threshold: #you MAY BE removing tumors with this sgID-bc
			if sgid_bc not in keepers:
				if contamination_report !="NA":
					contam_report_out.write("{},{},{},{},{},{},Rejected\n".format(sgid_bc, sgid, removal_threshold, pull_bc(sgid_bc), len(tumors_in[sgid_bc]), reads))
			else:
				for i in tumors_in[sgid_bc]: #you are now looping through occurences of the barcode
					sample = i.strip().split(",")[5]
					if sample == keepers[sgid_bc]: #this is the sample you're keeping the barcode in!
						outfile.write("{}\n".format(i))
						if contamination_report !="NA":
							contam_report_out.write("{},{},{},{},{},{},Rescued\n".format(sgid_bc, sgid, removal_threshold, pull_bc(sgid_bc), len(tumors_in[sgid_bc]), reads))

					else:
						if contamination_report !="NA":
							contam_report_out.write("{},{},{},{},{},{},Rejected\n".format(sgid_bc, sgid, removal_threshold, pull_bc(sgid_bc), len(tumors_in[sgid_bc]), reads))		
		else:
			for i in tumors_in[sgid_bc]:
				outfile.write("{}\n".format(i))
			if contamination_report !="NA":
				contam_report_out.write("{},{},{},{},{},{},Kept\n".format(sgid_bc, sgid, removal_threshold, pull_bc(sgid_bc), len(tumors_in[sgid_bc]), reads))

	if contamination_report != "NA":
		contam_report_out.close()
		



	#test = poisson.cdf(k=1,mu=_lambda)
	#print("test is {}. This is for n = {} and lambda = {}\n".format(test, 1, _lambda))

