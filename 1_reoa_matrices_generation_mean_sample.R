#############################
### REOA data preparation ###
#############################

### Load libraries and data
## libraries
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(ggplot2)
library(ggfortify)
library(limma)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(janitor)
library(dplyr)
library(tidyr)
library(mice)
library(reshape2)
library(openxlsx)
library(NbClust)
library(factoextra)
library(readr)
library(UniProt.ws)

## Data importation
## raw_data
raw_data <- read_csv("./raw_data/ITACAT_imp50_t.csv", col_names = TRUE)
# exp matrix
expression_mat <- as.data.frame(t(raw_data[-c(1:12)]))
colnames(expression_mat) <- raw_data$Accession
boxplot(expression_mat)

## Annotation
annotation <- openxlsx::read.xlsx("./raw_data/raw_data_ITACAT.xlsx")
annotation <- annotation[,c(1:3)]

### mean matrix generation
# mean
row_means <- rowMeans(expression_mat)
expression_mat <- cbind(expression_mat, mean_sample = row_means)

### Datasets for REOA
## expresson_matrix for REOA
# problem illness
index_problem_ill <- "C110"
index_problem_ill <- grepl(index_problem_ill, colnames(expression_mat))
expression_matrix_problem_illness <- expression_mat[,index_problem_ill,drop = F]
# write.table(x = expression_matrix_problem_illness,
# file = "./results/expression_mats/3_med_sample/expression_matrix_problem_illness.dat", sep = "\t",
#             col.names = F, row.names = F)
# problem control
index_problem_control <- "mean_sample" 
index_problem_control <- grepl(index_problem_control, colnames(expression_mat))
expression_matrix_problem_control <- expression_mat[,index_problem_control, drop = F] 
# write.table(x = expression_matrix_problem_control,
#  file = "./results/expression_mats/3_med_sample/expression_matrix_problem_control.dat", sep = "\t",
#              col.names = F, row.names = F)
# expression matrix for the
indexes_problem <- "C110|mean_sample"
indexes_problem <- grep(indexes_problem, colnames(expression_mat))
expression_matrix_cohort <- expression_mat[,-c(indexes_problem), drop = F]
# write.table(x = expression_matrix_cohort,
# file = "./results/expression_mats/3_med_sample/expression_matrix_to_reoa.dat", sep = "\t",
#             col.names = F, row.names = F)
 
### Annotation curation
## Accession_1 and Gene.names_1 generation
annotation <- annotation %>%
  rename(Protein.Groups = Accession,
         Protein.Description = Protein.names) %>%
  group_by(Protein.Groups) %>% 
  mutate(Accession_1 = unlist(strsplit(Protein.Groups, split = "\\;"),use.names = F)[1]) %>%
  mutate(Gene.names_1 = unlist(strsplit(Gene.names, split = "\\;"),use.names = F)[1]) %>% 
  ungroup() %>% 
  relocate(Accession_1, .after = Protein.Groups) %>% 
  relocate(Gene.names_1, .after = Gene.names)  

## Uniprot.ws data addition
# Uniprot.ws query
unip_query <- UniProt.ws::mapUniProt("UniProtKB_AC-ID", "UniProtKB", 
                                     query = annotation$Accession_1)
unip_query <- unip_query %>% 
  subset(select = c(From,Entry.Name,Protein.names,Gene.Names)) %>% 
  rename(
    Accession_1 = From,
    Protein.Name = Entry.Name,
    Protein.Description = Protein.names,
    Gene.names = Gene.Names) # %>% 
  # group_by(Accession_1) %>% 
  # mutate(Gene.Names_1 = unlist(strsplit(Gene.Names,split = " "),use.names = F)[1]) %>% 
  # ungroup()

# Uniprot.ws data addition
annotation <- annotation %>% 
  group_by(Protein.Groups) %>%
  # Add Protein.Name
  mutate(Protein.Name = ifelse(Accession_1 %in% unip_query$Accession_1, 
                unip_query$Protein.Name[match(Accession_1,
                                              unip_query$Accession_1)],NA)) %>% 
  # Add Protein.Description
  mutate(Protein.Description = ifelse(Accession_1 %in% unip_query$Accession_1, 
                                      unip_query$Protein.Description[match(Accession_1,
                                                                    unip_query$Accession_1)],NA)) %>% 
  # Add Gene.names if not in there
  # remove the space by ;
  # add Gene.names_1 if required
  mutate(Gene.names = ifelse(Accession_1 %in% unip_query$Accession_1 & is.na(Gene.names),
                                      unip_query$Gene.names[match(Accession_1,
                                                                           unip_query$Accession_1)],Gene.names)) %>%

  mutate(Gene.names = gsub(" ","\\;", Gene.names)) %>%
  mutate(Gene.names_1 = ifelse(is.na(Gene.names_1),
                               unlist(strsplit(Gene.names,split = "\\;"),use.names = F)[1],Gene.names_1)) %>% 
  
  relocate(Protein.Name, Protein.Description, .after = Gene.names_1)






