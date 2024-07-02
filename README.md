# reoa Downstream Analysis
## Credits
All the information in this README reffering to the execution of REOA can be found here:
https://github.com/pathint/reoa. All the credits must go to https://github.com/pathint the person/s who invented the approach and that have helped me and colegues to understand their application when it has been necessary.

## REOA cheatcode
1. Stable pairs calculation
    * reoa -j 0 1 599 -v expression_matrix_cohort.dat 6 1
    * mv stable_pairs_0.dat stable_pairs_0_verbose.dat 
    * reoa -j 0 1 599 expression_matrix_cohort.dat 6 1
2. Differentail expression analysis
    * bash onecomp -a 0 -m 0 -f 0.05 -v stable_pairs_0.dat expression_matrix_problem_control.dat expression_matrix_problem_illness.dat
    * bash cellcomp -a 0 -f 0.05 -v expression_matrix_problem_control.dat expression_matrix_problem_illness.dat
### Arguments
#### reoa, stable pairs calculation
-j 0: type of job to do the reoa, 0 for the stable-pairs calculations \
number of files for the stable_pairs calculation \
number of proteins/genes/CpGs to test \
-v: verbose results (add number of exceptions and pvalue) \
expression matrix for the cohort \
number of samples in the cohort sample \
if integer max number of exceptions for the calculations \
if float pvalue significance threshold for the stable pairs calculation

#### onecomp dysregulation state
-a: 0 Use RankCompv2, 1 Use RankCompv1 \
-m: 0 Use only one sample as problem_control, use the whole problem_control matrix as problem_control \
-f: significance threshold for the results \
-v: verbose results (add pvalue) \
stable_pairs result (NO verbose) \
expression_matrix_problem_control \
expression_matrix_problem_illness

#### cellcomp dysregulation state
TO DO !!!

## Vocabulary to be familiar with
TO DO !!!

## Rscripts for data preparation and data downstream analyiss
R code for the analysis of RankCompv2 data analysis for Onecomp applying for the analysis of a single sample against a Cohort of controls.
1. Scripts for different preprocessing versions
    * Starting with "1" data preprocessing unitil the REOA matricess are created and annotation generation. Run this script before runing script 2 and comment or remove the "write.table" function execution to avoid overwritting the raw_REOA_files.
    * Starting with "2" are for the downstream analysis of the output generated after running RankComp and Onecomp. 
2. Different  versions for the analysis and matrix generation
    * e72: A single sample is the "problem_control"
    * e72_c2_e32: Three samples are the "problem_control"
    * mean_sample: The mean sample of all the cohort samples is the "problem_control"
    * mean_sample_weight:  The weighted mean sample of all the cohort samples is the "problem_control"
    * most_similar_to_all: The sample with the highest accumulated index of correlation is the "problem_control"