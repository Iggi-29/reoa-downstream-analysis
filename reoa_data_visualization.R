###############################
### reoa data visualization ###
###############################

### Load libraries and data
# libraries
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
library(ggpattern)
library(ggnewscale)

#### Data importation
## raw_data - this will be an expression matrix
raw_data <- read_csv("./raw_data/ITACAT_imp50_t.csv", col_names = TRUE)

# exp matrix - remove unwanted columns
expression_mat <- as.data.frame(t(raw_data[-c(1:12)]))
colnames(expression_mat) <- raw_data$Accession
boxplot(expression_mat)

#### Recreation of the sets that REOA has processed 
# Extract the data that has already been processed by reoa, 
# in the future source another R script to do this 
## problem - ill sample
expression_matrix_C110 <- expression_mat[,c("C110"), drop = F]
## problem - control sample
expression_matrix_E72 <- expression_mat[,c("E72"), drop = F] ### THIS IS THE SAMPLE OF S93 first collection !!
## cohort
expression_matrix_reoa <- expression_mat[,-c(23,51), drop = F]

#### Annotation
annotation <- openxlsx::read.xlsx("./raw_data/raw_data_ITACAT.xlsx")
annotation <- annotation[,c(1:3)]

### Uniprot.ws work
annotation <- annotation %>%
  group_by(Accession) %>% 
  mutate(Accession_1 = unlist(strsplit(Accession, split = "\\;"),use.names = F)[1]) %>% 
  ungroup() %>% 
  relocate(Accession_1, .after = Accession)

annotation_uniP <- UniProt.ws::mapUniProt("UniProtKB_AC-ID", "UniProtKB", 
                                          query = annotation$Accession_1)
annotation <- annotation %>% 
  group_by(Accession) %>% 
  mutate(Gene.names = ifelse(is.na(Gene.names),
                             annotation_uniP$Gene.Names[match(Accession_1, annotation_uniP$From)],
                             Gene.names)) %>% 
  mutate(Gene.names = gsub(x =Gene.names, pattern = " ", replacement = ";")) %>% 
  subset(select = -c(Accession_1))

#### Import REOA data
## !! ALWAYS SUM+1 as in reoa 0 means "first row" and in R, this information would be 1
## Stable pairs
REOs <- readr::read_tsv("./raw_data/reoa_pval_001_c110_e72/stable_pairs_0.dat", col_names = F)
colnames(REOs) <- c("high","low","Exeptions","Pvalue")
REOs$high <- REOs$high+1
REOs$low <- REOs$low+1
REOs_genes <- unique(c(REOs$high,REOs$low))

## Upregulated
REOs_up <- readr::read_tsv("./raw_data/reoa_pval_001_c110_e72/up_regulated_0.dat", col_names = F)
colnames(REOs_up) <- c("Protein","FDR")
REOs_up$Protein <- REOs_up$Protein+1

## Downregulated
REOs_down <- readr::read_tsv("./raw_data/reoa_pval_001_c110_e72/down_regulated_0.dat", col_names = F)
colnames(REOs_down) <- c("Protein","FDR")
REOs_down$Protein <- REOs_down$Protein+1

#### Subset the expression matrices to gather the data for problem-ill/problem-control and cohort
## Cohort expression matrix
expression_matrix <- expression_mat
expression_matrix$Accession <- row.names(expression_matrix)
expression_matrix <- merge(annotation, expression_matrix, by = "Accession")

individualized_cohort <- expression_matrix
individualized_cohort <- individualized_cohort %>% 
  subset(select = -c(C110,E72))
samples <- colnames(individualized_cohort)[c(4:ncol(individualized_cohort))]

# To long format
individualized_cohort <- individualized_cohort %>% 
  tidyr::pivot_longer(cols = all_of(samples), 
                      names_to = "sample", values_to = "abundance")

## Problem control sample 
### Get the matrix
individualized_control_samp <- expression_matrix_E72
individualized_control_samp$Accession <- rownames(individualized_control_samp)
rownames(individualized_control_samp) <- NULL
individualized_control_samp$number <- rownames(individualized_control_samp)
individualized_control_samp <- merge(annotation, individualized_control_samp,
                                     by = "Accession")
### mimic the cohort data
individualized_control_samp <- individualized_control_samp %>% 
  mutate(condition = ifelse(number %in% REOs_up$Protein,"Upregulated",
                            ifelse(number %in% REOs_down$Protein,"Downregulated","Non-Dysregulated"))) %>% 
  mutate(condition = ifelse(!number %in% REOs_genes,
                            "Non-stablepair",condition))

individualized_control_samp$sample <- "E72"
colnames(individualized_control_samp)[4] <- "abundance"
individualized_control_samp <- individualized_control_samp %>% 
  relocate(number, .after = condition)

## Problem ill sample
# Get the matrix
individualized_ill_samp <- expression_matrix_C110
individualized_ill_samp$Accession <- rownames(individualized_ill_samp)
rownames(individualized_ill_samp) <- NULL
individualized_ill_samp$number <- rownames(individualized_ill_samp)
individualized_ill_samp <- merge(annotation, individualized_ill_samp,
                                 by = "Accession")
# mimic the cohort data
individualized_ill_samp <- individualized_ill_samp %>% 
  mutate(condition = ifelse(number %in% REOs_up$Protein,"Upregulated",
                            ifelse(number %in% REOs_down$Protein,"Downregulated","Non-Dysregulated"))) %>% 
  mutate(condition = ifelse(!number %in% REOs_genes,
                            "Non-stablepair",condition))

individualized_ill_samp$sample <- "C110"
colnames(individualized_ill_samp)[4] <- "abundance"
individualized_ill_samp <- individualized_ill_samp %>% 
  relocate(number, .after = condition)

# Cohort matrix w/up and down info
individualized_cohort <- individualized_cohort %>% 
  mutate(condition = ifelse(Gene.names %in% 
                              individualized_ill_samp$Gene.names,
                            individualized_ill_samp$condition[match(Gene.names,individualized_ill_samp$Gene.names)],
                            NA)) %>% 
  mutate(number = ifelse(Gene.names %in% 
                           individualized_ill_samp$Gene.names,
                         individualized_ill_samp$number[match(Gene.names,individualized_ill_samp$Gene.names)],
                         NA))

### iDEA with unified info for the 3 samples
iDEA <- rbind(individualized_cohort, individualized_control_samp,individualized_ill_samp)

#### Downstream Analysis for the IDEA
## Create a meta_sata
meta_data <- as.data.frame(tibble::tibble(
  sample = c(samples,"E72","C110"),
  type_of_sample = c(rep("Control Cohort", times = length(samples)),
                     "Problem Control","Problem Ill")
))

# Add meta_data info to the iDEA dataframe
iDEA <- iDEA %>% 
  mutate(sample_type = ifelse(sample %in% meta_data$sample,
                              meta_data$type_of_sample[match(sample, 
                                                             meta_data$sample)],NA))
## Mutate Gene.nanmes
iDEA <- iDEA %>% 
  group_by(sample) %>% 
  filter(!is.na(Gene.names)) %>% 
  filter(!duplicated(Gene.names)) %>%
  ungroup() %>% 
  group_by(Accession) %>% 
  mutate(Gene.names_1 = unlist(strsplit(Gene.names, split = "\\;"),use.names = F)[1]) %>% 
  ungroup() %>% 
  relocate(Gene.names_1, .after = Gene.names)

## Add the pvalue info for the REOA analysis 
REOs_upn_down <- rbind(REOs_up, REOs_down)
REOs_upn_down$Protein <- as.character(REOs_upn_down$Protein)
iDEA <- iDEA %>% 
  mutate(REOA_pval = ifelse(number %in% REOs_upn_down$Protein,
                            REOs_upn_down$FDR[match(number,
                                                    REOs_upn_down$Protein)],1))
### Finish the iDEA dataset
# Calculate the median
iDEA_to_plot <- iDEA %>%
  group_by(sample_type,Gene.names) %>% 
  mutate(MED = mean(abundance ,na.rm = T)) %>% 
  ungroup() %>% 
  subset(select = -c(sample,abundance)) %>% 
  distinct()

# Calculate RANK
iDEA_to_plot <- iDEA_to_plot %>% 
  group_by(sample_type) %>% 
  mutate(rank = rank(MED)) %>% 
  ungroup()

# Remove redundant information
iDEA_genes <- iDEA_to_plot %>% 
  subset(select = c(Gene.names, REOA_pval, condition)) %>% 
  distinct() %>% 
  mutate(condition = factor(condition, levels = c("Downregulated","Upregulated","Non-Dysregulated"))) %>%
  arrange(condition,REOA_pval)
genes <- iDEA_genes$Gene.names

#### Information for the plots
# Set the order
iDEA_to_plot <- iDEA_to_plot %>%
  mutate(condition = factor(condition, levels = c("Downregulated","Non-stablepair","Non-Dysregulated","Upregulated"))) %>%
  mutate(sample_type = factor(sample_type, levels = c("Control Cohort","Problem Control","Problem Ill"))) %>% 
  mutate(Gene.names = factor(Gene.names, levels = c(genes))) %>% 
  arrange(condition,sample_type,Gene.names)

# Create a new combined factor column to ensure the good ordering of the data
iDEA_to_plot$combined <- factor(paste(iDEA_to_plot$condition, 
                                      iDEA_to_plot$Gene.names, 
                                      iDEA_to_plot$sample_type, sep = "_"), 
                                levels = unique(paste(iDEA_to_plot$condition, 
                                                      iDEA_to_plot$Gene.names, 
                                                      iDEA_to_plot$sample_type, sep = "_")))
# Set colors for the graphs
types_of_samples <- c("Control Cohort" = "gray",
                      "Problem Control" = "purple",
                      "Problem Ill" = "yellow")
dynamics_ <- c("Downregulated" = "red",
               "Upregulated" = "blue",
               "Non-Dysregulated" = "black",
               "Non-stablepair" = "darkgreen")

#### Barplots
### General Barplot
REOA_bars <- ggplot(iDEA_to_plot, aes(y=MED, x=combined))+
  scale_x_discrete(limits = levels(iDEA_to_plot$combined),
                   labels = iDEA_to_plot$Gene.names)+
  geom_bar(mapping = aes(color = condition), 
           position="dodge", stat="identity")+
  scale_color_manual(values = dynamics_)+
  theme_minimal() +
  labs(title = "REOA results",
       x = "Protein")+
  theme_minimal()+
  theme(axis.title = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())+
  guides(size = "none")
REOA_bars

### Do a barplot for the top 5 up and downregulated proteins 
# Get the genes
up_to_check <- REOs_up$Protein[order(REOs_up$FDR)[1:5]]
down_to_check <- REOs_down$Protein[order(REOs_down$FDR)[1:5]]

## Get the data
iDEA_to_plot2 <- iDEA_to_plot %>% 
  dplyr::filter(number %in% c(up_to_check,down_to_check))

iDEA_to_plot2 <- iDEA_to_plot2 %>%
  mutate(condition = as.character(condition)) %>% 
  mutate(sample_type = as.character(sample_type)) %>% 
  mutate(Gene.names = as.character(Gene.names))  

iDEA_to_plot2 <- iDEA_to_plot2 %>%
  mutate(condition = factor(condition, levels = c("Downregulated","Non-stablepair","Non-Dysregulated","Upregulated"))) %>%
  mutate(sample_type = factor(sample_type, levels = c("Control Cohort","Problem Control","Problem Ill"))) %>% 
  mutate(Gene.names = factor(Gene.names, levels = c(genes))) %>% 
  arrange(condition,Gene.names,sample_type)

iDEA_to_plot2$combined <- factor(paste(iDEA_to_plot2$condition, 
                                       iDEA_to_plot2$Gene.names, 
                                       iDEA_to_plot2$sample_type, sep = "_"), 
                                 levels = unique(paste(iDEA_to_plot2$condition, 
                                                       iDEA_to_plot2$Gene.names, 
                                                       iDEA_to_plot2$sample_type, sep = "_")))
iDEA_to_plot2 <- iDEA_to_plot2 %>% 
  mutate(Stripe = ifelse(condition == "Non-stablepair","stripe","none"))

# old plot Model
REOA_bars2 <- ggplot(iDEA_to_plot2, aes(y=MED, x=combined))+
  scale_x_discrete(limits = levels(iDEA_to_plot2$combined),
                   labels = iDEA_to_plot2$Gene.names)+
  geom_bar(mapping = aes(fill = sample_type, color = condition), 
           position="dodge", stat="identity")+
  scale_fill_manual(values = types_of_samples)+
  scale_color_manual(values = dynamics_)+
  theme_minimal() +
  labs(title = "REOA results",
       x = "Protein")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(size = "none")

## Barplot
REOA_bars2 <- ggplot(iDEA_to_plot2, aes(y=MED, x=combined))+
  scale_x_discrete(limits = levels(iDEA_to_plot2$combined),
                   labels = iDEA_to_plot2$Gene.names)+
  geom_bar(mapping = aes(fill = sample_type, color = condition), 
           position="dodge", stat="identity", linewidth = 1)+
  scale_fill_manual(values = types_of_samples)+
  scale_color_manual(values = dynamics_)+
  geom_bar_pattern(fill = "white",alpha = 0, aes(pattern = Stripe, color = condition),
                   position = "dodge", stat = "identity",linewidth=1)+
  scale_pattern_identity()+
  theme_minimal() +
  labs(title = "REOA results\nTop5 Up and Down proteins",
       x = "Protein")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(size = "none")

REOA_bars2


### Do a barplot to an upregulated pathway  
# Get the genes
up_to_check <- openxlsx::read.xlsx("./results/BP_C110_E72_UP_pval001.xlsx")
up_to_check <- up_to_check %>% 
  arrange(desc(Count)) %>% 
  slice_head(n = 1)
up_to_check <- up_to_check$geneID
up_to_check <- unlist(strsplit(up_to_check, split = "\\/"), use.names = F)

## Get the data
iDEA_to_plot3 <- iDEA_to_plot %>% 
  dplyr::filter(Gene.names_1 %in% c(up_to_check))

iDEA_to_plot3 <- iDEA_to_plot3 %>%
  mutate(condition = as.character(condition)) %>% 
  mutate(sample_type = as.character(sample_type)) %>% 
  mutate(Gene.names = as.character(Gene.names))  

iDEA_to_plot3 <- iDEA_to_plot3 %>%
  mutate(condition = factor(condition, levels = c("Downregulated","Non-stablepair","Non-Dysregulated","Upregulated"))) %>%
  mutate(sample_type = factor(sample_type, levels = c("Control Cohort","Problem Control","Problem Ill"))) %>% 
  mutate(Gene.names = factor(Gene.names, levels = c(genes))) %>% 
  arrange(condition,sample_type,Gene.names)

iDEA_to_plot3$combined <- factor(paste(iDEA_to_plot3$condition, 
                                       iDEA_to_plot3$Gene.names, 
                                       iDEA_to_plot3$sample_type, sep = "_"), 
                                 levels = unique(paste(iDEA_to_plot3$condition, 
                                                       iDEA_to_plot3$Gene.names, 
                                                       iDEA_to_plot3$sample_type, sep = "_")))

iDEA_to_plot3 <- iDEA_to_plot3 %>% 
  mutate(Stripe = ifelse(condition == "Non-stablepair","stripe","none"))

## Barplot
REOA_bars3 <- ggplot(iDEA_to_plot3, aes(y=MED, x=combined))+
  scale_x_discrete(limits = levels(iDEA_to_plot3$combined),
                   labels = iDEA_to_plot3$Gene.names)+
  geom_bar(mapping = aes(fill = sample_type, color = condition), 
           position="dodge", stat="identity", linewidth = 1)+
  scale_fill_manual(values = types_of_samples)+
  scale_color_manual(values = dynamics_)+
  geom_bar_pattern(fill = "white",alpha = 0, aes(pattern = Stripe, color = condition),
                   position = "dodge", stat = "identity",linewidth=1)+
  scale_pattern_identity()+
  theme_minimal() +
  labs(title = "REOA results\nTo an upregulated pathway",
       x = "Protein")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(size = "none")

REOA_bars3


### Do a barplot to a randomly selected pathway
# Get the genes
up_to_check <- openxlsx::read.xlsx("./results/BP_C110_E72_UP.xlsx")
up_to_check <- up_to_check %>%
  filter(Description == "complement activation, classical pathway")
up_to_check <- up_to_check$geneID
up_to_check <- unlist(strsplit(up_to_check, split = "\\/"), use.names = F)

## Get the data
iDEA_to_plot4 <- iDEA_to_plot %>% 
  dplyr::filter(Gene.names_1 %in% c(up_to_check))

iDEA_to_plot4 <- iDEA_to_plot4 %>%
  mutate(condition = as.character(condition)) %>% 
  mutate(sample_type = as.character(sample_type)) %>% 
  mutate(Gene.names = as.character(Gene.names))  

iDEA_to_plot4 <- iDEA_to_plot4 %>%
  mutate(condition = factor(condition, levels = c("Downregulated","Non-stablepair","Non-Dysregulated","Upregulated"))) %>%
  mutate(sample_type = factor(sample_type, levels = c("Control Cohort","Problem Control","Problem Ill"))) %>% 
  mutate(Gene.names = factor(Gene.names, levels = c(genes))) %>% 
  arrange(condition,sample_type,Gene.names)

iDEA_to_plot4$combined <- factor(paste(iDEA_to_plot4$condition, 
                                       iDEA_to_plot4$Gene.names, 
                                       iDEA_to_plot4$sample_type, sep = "_"), 
                                 levels = unique(paste(iDEA_to_plot4$condition, 
                                                       iDEA_to_plot4$Gene.names, 
                                                       iDEA_to_plot4$sample_type, sep = "_")))

iDEA_to_plot4 <- iDEA_to_plot4 %>% 
  mutate(Stripe = ifelse(condition == "Non-stablepair","stripe","none"))

## Barplot
REOA_bars4 <- ggplot(iDEA_to_plot4, aes(y=MED, x=combined))+
  scale_x_discrete(limits = levels(iDEA_to_plot4$combined),
                   labels = iDEA_to_plot4$Gene.names)+
    geom_bar(mapping = aes(fill = sample_type, color = condition), 
             position="dodge", stat="identity", linewidth = 1)+
  scale_fill_manual(values = types_of_samples)+
  scale_color_manual(values = dynamics_)+
  geom_bar_pattern(fill = "white",alpha = 0, aes(pattern = Stripe, color = condition),
                   position = "dodge", stat = "identity",linewidth=1)+
  scale_pattern_identity()+
  theme_minimal() +
  labs(title = "REOA results\nWith a 'randomly' selected pathway",
       x = "Protein")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(size = "none")
  
REOA_bars4

#### Volcano plot
## Get the data
iDEA_to_plot5 <- iDEA_to_plot 

# Remove 'Problem Control' samples and Non-stable and Non-Dysregulated proteins
iDEA_to_plot5 <- iDEA_to_plot5 %>%
  mutate(condition = as.character(condition)) %>% 
  mutate(sample_type = as.character(sample_type)) %>% 
  mutate(Gene.names = as.character(Gene.names)) %>% 
  filter(sample_type != "Problem Control") %>% 
  filter(!condition %in% c("Non-stablepair","Non-Dysregulated")) 

# Get the Fold change between the cohort and the 'Problem Ill' protein
wide_idea_to_plot5 <- iDEA_to_plot5 %>% 
  tidyr::spread(sample_type, MED)

wide_idea_to_plot5 <- wide_idea_to_plot5 %>% 
  subset(select = c(Gene.names,`Control Cohort`,`Problem Ill`)) %>% 
  group_by(Gene.names) %>% 
  mutate(`Control Cohort` = min(`Control Cohort`, na.rm = T)) %>% 
  mutate(`Problem Ill` = min(`Problem Ill`, na.rm = T)) %>% 
  ungroup() %>% 
  distinct() %>% 
  mutate(logFC = `Problem Ill` - `Control Cohort`)

iDEA_to_plot5 <- iDEA_to_plot5 %>% 
  mutate(Ill_vs_Cohort = wide_idea_to_plot5$logFC[match(Gene.names,
                                                wide_idea_to_plot5$Gene.names)]) %>% 
  subset(select = c(Accession,Protein.names,Gene.names,Gene.names_1,REOA_pval,Ill_vs_Cohort)) %>% 
  distinct()

# The 0 PValues
pval_no_0 <- sort(unique(iDEA_to_plot5$REOA_pval))[2] 

iDEA_to_plot5 <- iDEA_to_plot5 %>% 
  mutate(REOA_pval = ifelse(REOA_pval == 0,pval_no_0,
                            REOA_pval))

### Volcano plot
# data manipulation
raw_data <- raw_data %>% 
  filter(!is.na(Gene.names)) %>% 
  filter(!is.na(Protein.names))

raw_data$fold_change <- ifelse(test = raw_data$logFC > 0, 
                               yes = 2^raw_data$logFC,
                               no = (-1/2^raw_data$logFC))
raw_data <- relocate(.data = raw_data, fold_change, .after = logFC)
raw_data$fold_change <- raw_data$fold_change*-1
raw_data$logFC <- raw_data$logFC*-1

raw_data <- raw_data %>% 
  mutate(gene_type = ifelse(fold_change >= 0.3785116 & adj.P.Val <= 0.05, "Enriched in G1",
                            ifelse(fold_change <= -0.3785116 & adj.P.Val <= 0.05,"Enriched in G2","non-enriched")))

# Define some variables
cols <- c("Enriched in G2" = "#ff00ff", "Enriched in G1" = "#ffff00", "non-enriched" = "grey")


# plotting
vl_plot <- ggplot(data = raw_data,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene.names))+ 
  geom_point(aes(x = logFC, y = -log10(adj.P.Val),
                 colour = (gene_type),
                 size = GRAN))+
  scale_colour_manual(values = cols)+
  scale_size_manual(values = c(6,10))+
  geom_point(data = raw_data[raw_data$GRAN == "SI",],
             pch=21, fill=NA, size=10, colour="black", stroke = 2.5)+
  
  geom_text_repel(data = . %>% 
                    filter(raw_data$GRAN %in% c("SI")) %>%
                    filter(logFC > 0) %>% 
                    filter(!Gene.names %in% c("HBB","RDH11","HBA1","P2RX1","STXBP3","GNAQ",
                                              "UBASH3B","GSTO1","ILK","TXN","PRDX1","GNA13")),
                  nudge_x = 0.50, 
                  nudge_y = .1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 2)+
  
  geom_text_repel(data = . %>% 
                    filter(raw_data$GRAN %in% c("SI")) %>%
                    filter(logFC > 0) %>% 
                    filter(Gene.names %in% c("HBB","RDH11","HBA1","P2RX1","STXBP3","GNAQ",
                                             "UBASH3B","GSTO1","ILK","TXN","PRDX1","GNA13")),
                  nudge_x = -0.50, 
                  nudge_y = .1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 2)+
  
  geom_text_repel(data = . %>% 
                    filter(raw_data$GRAN %in% c("SI")) %>% 
                    filter(logFC < 0), 
                  nudge_x = -1, 
                  nudge_y = 0,
                  segment.ncp = 3,
                  segment.angle = 2)+
  
  ####  geom_text_repel(data = . %>% 
  ####                    filter(raw_data$GRAN %in% c("SI")) %>%
  ####                    filter(logFC > 0) %>%
  ####                    filter(Gene.names %in% c("HBB","TXNDC17","UBE2V2"), .preserve = T),
  ####                  nudge_x = 1, 
  ####                  nudge_y = 0.2,
  ####                  segment.ncp = 3,
  ####                  segment.angle = 2)+
  
  
geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(0.3785116, -0.3785116), linetype = "dashed")+
  scale_x_continuous(breaks = seq(-4,4,2),
                     limits = c(-4,4))+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))+
  xlab("G1 vs. G2 (LogFC)")
vl_plot





