# reoa Downstream Analysis
## REOA cheatcode
1. Stable pairs calculation
    * reoa -j 0 1 599 -v expression_matrix_cohort.dat 6 1
    * mv stable_pairs_0.dat stable_pairs_0_verbose.dat 
    * reoa -j 0 1 599 expression_matrix_cohort.dat 6 1
2. Differentail expression analysis
    * bash onecomp -a 0 -m 0 -f 0.05 -v stable_pairs_0.dat expression_matrix_problem_control.dat expression_matrix_problem_illness.dat
    * bash cellcomp -a 0 -f 0.05 -v expression_matrix_problem_control.dat expression_matrix_problem_illness.dat

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