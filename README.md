# Aging
Quantifying the impact of aging on lung tumorigenesis and tumor suppressor inactivation

This repository contains code associated with the manuscript

**Shuldiner E. G.**, **Karmakar S.**, **Tsai M. K.**, **Hebert J. D.**, **Tang Y. J.**, **Andrejka L.**, **Robertson M. R.**, **Wang M.**, **Detrick C. R.**, **Cai H.**, **Tang R.**, **Kunder C. A.**, **Feldser D. M.**, **Petrov D. A.**, and **Winslow M. M.** _Aging represses oncogenic KRAS-driven lung tumorigenesis and alters tumor suppression_. Currently available at https://www.biorxiv.org/content/10.1101/2024.05.28.596319v1

# Running the tuba-seq pipeline

## Environment
Dependencies required to run the code can be found in ./Environment. Create the virtual environment and load it using conda.

```
conda create -f Environment/tubaseq.yml
source activate tubaseq
```

## Project file

Information defining a project (including sample names, locations, and characteristics, as well as the makeup of the viral pool) are passed through a "project file". The project file is tab-separated, with each line of the file defining a parameter. Parameter values with multiple fields are comma-separated; see example project file.

By default, the location and name of the project file should be as follows:
<root>/Parameters/<project_id>_project_file.txt

The following parameters must be defined in the project file (if a parameter is not relevant to the analysis, substitute with NA):

| Parameter Name    | Definition | Datatype |
| -------- | ------- | ------- |
| PROJECT  | Name of the project  | string |
| SAMPLES | Names of samples included in the project    | comma-separated strings |
| SAMPLES_PATHS    | Full paths to R1 fastq files for; R2 must be in the same directory. Indexing should match <SAMPLES>.   | comma-separated strings |
| GENOTYPES | Genotypes of each sample; indexing should match <SAMPLES>. | comma-separated strings |
| TREATMENT | Treatments of each sample (e.g. drug treatment, age, etc); indexing should match <SAMPLES> | comma-separated strings |
| TITERS | Titer of virus delivered to each sample; indexing should match <SAMPLES> | comma-separated integers |
| SPIKE_BC | Barcodes of spike-in cell lines associated with each sample | Barcodes are colon-separated; values for each sample should be either comma-separated with indexing matching <SAMPLES> or provide a single value to be assumed for all samples |
|SIZE_SPIKE | Number of cells for each spike-in cell line added to each sample | Values for each barcode are colon separated and indexing should match <SPIKE_BC>. Values for each sample should be comma separated with indexing matching <SAMPLES> or provide a single value to be assumed for all samples. |
| SGIDS | sgRNA-identifying sequences | comma-separated strings |
| SGRNAS | Names/identifiers for sgRNAs used in the study; indexing should match <SGIDS>. | comma-separated strings |
| GENES | Names of genes targeted by sgRNAs used in the study; indexing should match <SGIDS>. | 
| INERT | Non- or safe-targeting control sgRNAs (used as baseline in calculating the effects of gene inactivation) | comma-separated strings |
| CAS9NEG_GT | Genotype of mice lacking Cas9 (matching one of <GENOTYPES>) | string |

## Parameter file
The parameter values used in performing a given analysis are passed through a "parameter file". The parameter file is tab-separated, with each line of the file defining a parameter.

By default, the location and name of the project file should be as follows:
<root>/Parameters/<parameter_id>_parameter_file.txt

The following parameters must be defined in the project file (see Methods section of manuscript for complete definition):

| Parameter Name    | Definition | Datatype |
| -------- | ------- | ------- |
| PARAMETER_SET  | Identifier for the set of parameters  | string |
| DIST_BETWEEN_BC_FOR_READ | Allowable hamming distance between barcode sequence on R1 and R2   | integer |
| R1_REGEX_LOOSENESS    | Allowable error in matching of fuzzy regex for barcode on R1.   | integer |
| R2_REGEX_LOOSENESS | Allowable error in matching of fuzzy regex for barcode on R2 | integer |
| BC_SIMILARITY_THRESHOLD | Similarity (hamming distance) threshold for identifying 'spurious tumors' likely to have arise due to PCR or sequencing error. | integer |
| READ_ERROR_SIZE_THRESHOLD | Percent of size of larger tumor to collapse smaller tumor if barcodes have a hamming distance within <BC_SIMILARITY_THRESHOLD>| float |
| CONTAMINATION_REMOVAL_THRESHOLD | Tail probability used to define threshold Nr for identification of instances of barcode recurrence unlikely to be caused by chance. | float |

An example project and parameter file are provided in /Parameters. Note that pathes in the project file will need to be substituted. 

# Adaptive Sampling
In comparing the effects of tumor suppressor inactivation in young and old mice, we sought to compare equivalent portions of the tumor size distributions (i.e., the same number of tumors per infectious unit of virus delivered). To ensure this, we scaled the number of tumors analyzed for each sgRNA in each cohort to account for differences in viral titer and the number of mice transduced in each cohort, and then analyzed the largest tumors per Lenti-sgRNA/Cre vector. The manuscript contains additional details on these statistics and their interpretation; for code related to the implementation of adaptive sampling, see https://github.com/eshuldiner/Eml4-Alk.
