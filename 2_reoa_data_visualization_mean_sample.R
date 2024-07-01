##########################
### REOA Data analysis ###

#### Load libraries and data
# libraries
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(ggplot2)
library(ggfortify)
library(limma)
library(gplots)
# library()
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
library(scales)
library(clusterProfiler)
library(ggrepel)
library(ComplexHeatmap)

### Data importation
# import the matrices sourcing the code for the expression matrices and the annotation sets
source("./1_reoa_matrices_generation_mean_sample.R")
source("./functions/collapseGO.R")
source("./functions/minestrone.R")
source("./functions/treemaping.R")

###############
## IMPORTANT ##
###############
# Try to work with protein groups as unit Protein.Groups

#### Import REOA data 
## !! ALWAYS SUM+1 as in reoa 0 means "first row" and in R, this information would be 1
## Stable pairs
REOs <- readr::read_tsv("./results/reoa/3_med_sample/stable_pairs_0_verbose.dat", col_names = F)
colnames(REOs) <- c("high","low","Exeptions","Pvalue")
REOs$high <- REOs$high+1
REOs$low <- REOs$low+1
REOs_genes <- unique(c(REOs$high,REOs$low))
REOs$high <- as.character(REOs$high)
REOs$low <- as.character(REOs$low)

## Upregulated
REOs_up <- readr::read_tsv("./results/reoa/3_med_sample/up_regulated_0.dat", col_names = F)
colnames(REOs_up) <- c("Protein","FDR")
REOs_up$Protein <- REOs_up$Protein + 1
REOs_up$Protein <- as.character(REOs_up$Protein)

## Downregulated
REOs_down <- readr::read_tsv("./results/reoa/3_med_sample/down_regulated_0.dat", col_names = F)
colnames(REOs_down) <- c("Protein","FDR")
REOs_down$Protein <- REOs_down$Protein+1
REOs_down$Protein <- as.character(REOs_down$Protein)

# bind up and downregulated proteins info
REOs_upn_down <- rbind(REOs_up, REOs_down)

#### Separate information of the cohort and the problem_control/ill samples
## Cohort expression matrix
# samples
samples <- colnames(expression_matrix_cohort)

# work the expression matrix
expression_matrix_cohort$Protein.Groups <- rownames(expression_matrix_cohort)
rownames(expression_matrix_cohort) <- NULL
expression_matrix_cohort$number <- rownames(expression_matrix_cohort)
expression_matrix_cohort <- merge(annotation, expression_matrix_cohort, by = "Protein.Groups")

# To long format
individualized_cohort <- expression_matrix_cohort %>% 
  tidyr::pivot_longer(cols = all_of(samples), 
                      names_to = "sample", values_to = "abundance")
individualized_cohort$type_of_samples <- "Cohort" 
individualized_cohort <- individualized_cohort %>%
  group_by(Gene.names) %>%
  mutate(MED = mean(abundance)) %>%
  ungroup()

## Problem Control sample/s
samples <- colnames(expression_matrix_problem_control)

expression_matrix_problem_control$Protein.Groups <- rownames(expression_matrix_problem_control)
rownames(expression_matrix_problem_control) <- NULL
expression_matrix_problem_control$number <- rownames(expression_matrix_problem_control)
expression_matrix_problem_control <- merge(annotation, expression_matrix_problem_control, by = "Protein.Groups")

# To long format
individualized_problem_control <- expression_matrix_problem_control %>% 
  tidyr::pivot_longer(cols = all_of(samples), 
                      names_to = "sample", values_to = "abundance")
individualized_problem_control$type_of_samples <- "Problem Control" 
individualized_problem_control <- individualized_problem_control %>%
  group_by(Gene.names) %>%
  mutate(MED = abundance) %>%
  ungroup()

## Problem Ill smaple/s
samples <- colnames(expression_matrix_problem_illness)

expression_matrix_problem_illness$Protein.Groups <- rownames(expression_matrix_problem_illness)
rownames(expression_matrix_problem_illness) <- NULL
expression_matrix_problem_illness$number <- rownames(expression_matrix_problem_illness)
expression_matrix_problem_illness <- merge(annotation, expression_matrix_problem_illness, by = "Protein.Groups")

# To long format
individualized_problem_illness <- expression_matrix_problem_illness %>% 
  tidyr::pivot_longer(cols = all_of(samples), 
                      names_to = "sample", values_to = "abundance")
individualized_problem_illness$type_of_samples <- "Problem Illness" 
individualized_problem_illness <- individualized_problem_illness %>%
  group_by(Gene.names) %>%
  mutate(MED = abundance) %>%
  ungroup()

#### iDEA 
### Generation of iDEA dataframe with unified info for the 3 samples
iDEA <- rbind(individualized_cohort, individualized_problem_control,individualized_problem_illness)

## Add the up/downregulation iformation
iDEA <- iDEA %>%
  # Add the condition info 
  mutate(condition = ifelse(number %in% REOs_up$Protein,"Upregulated",
                            ifelse(number %in% REOs_down$Protein,"Downregulated","Non-Dysregulated"))) %>% 
  mutate(condition = ifelse(!number %in% REOs_genes,
                            "Non-stablepair",condition)) %>% 
  # Add the REOA_pvalue info
  mutate(REOA_pval = ifelse(number %in% REOs_upn_down$Protein,
                            REOs_upn_down$FDR[match(number,
                                                    REOs_upn_down$Protein)],1)) %>% 
  # If required, remove Empty, NA and duplicated Gene.names
  group_by(type_of_samples) %>% 
  filter(Gene.names != "") %>% 
  filter(Gene.names != " ") %>% 
  filter(!is.na(Gene.names)) %>%
  filter(!duplicated(Gene.names)) %>% 
  ungroup()
  
## Calculations for the REOA information
iDEA_to_plot <- iDEA %>%
  # For the Cohort samples, calculate the Median abundance (MED)
  # group_by(type_of_samples,Gene.names) %>%
  # mutate(MED = ifelse(type_of_samples == "Cohort",
  #                     mean(abundance),abundance)) %>%
  # mutate(MED = mean(abundance, rm.na = T)) %>%
  # ungroup() %>%
  # Remove the sample and abundance columns and then remove the redundant columns
  subset(select = -c(sample,abundance)) %>% 
  distinct() %>% 
  # Calculate RANK of the expression per sample
  group_by(type_of_samples) %>% 
  mutate(rank = rank(MED)) %>% 
  ungroup()

# Save the genes column (each gene only once) and in a arranged way
iDEA_genes <- iDEA_to_plot %>% 
  subset(select = c(Gene.names, REOA_pval, condition)) %>% 
  distinct() %>% 
  mutate(condition = factor(condition, levels = c("Downregulated","Upregulated","Non-Dysregulated"))) %>%
  arrange(condition,REOA_pval)
iDEA_genes <- iDEA_genes$Gene.names

#### Plotting
### General barplot
## Prepare the data
# Set the order of the entries of the dataframe in the plot
iDEA_to_plot <- iDEA_to_plot %>%
  mutate(condition = factor(condition, levels = c("Downregulated","Non-stablepair","Non-Dysregulated","Upregulated"))) %>%
  mutate(type_of_samples = factor(type_of_samples, levels = c("Cohort","Problem Control","Problem Illness"))) %>% 
  mutate(Gene.names = factor(Gene.names, levels = c(iDEA_genes))) %>% 
  arrange(condition,type_of_samples,Gene.names)

# Create a new combined factor column to ensure the good ordering of the data
iDEA_to_plot$combined <- factor(paste(iDEA_to_plot$condition, 
                                      iDEA_to_plot$Gene.names, 
                                      iDEA_to_plot$type_of_samples, sep = "_"), 
                                levels = unique(paste(iDEA_to_plot$condition, 
                                                      iDEA_to_plot$Gene.names, 
                                                      iDEA_to_plot$type_of_samples, sep = "_")))
# Set colors for the graphs
type_of_samples <- c("Cohort" = "gray",
                      "Problem Control" = "purple",
                      "Problem Illness" = "yellow")
dynamics_ <- c("Downregulated" = "red",
               "Upregulated" = "blue",
               "Non-Dysregulated" = "black",
               "Non-stablepair" = "darkgreen")

## Do the plot
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
ggsave(REOA_bars, filename = "./plots/3_med_sample/REOA_barmap.pdf")

### Do a barplot for the top 5 up and downregulated proteins 
# Get the genes
up_to_check <- REOs_up$Protein[order(REOs_up$FDR)[1:5]]
down_to_check <- REOs_down$Protein[order(REOs_down$FDR)[1:5]]

## Get the data
iDEA_to_plot2 <- iDEA_to_plot %>% 
  dplyr::filter(number %in% c(up_to_check,down_to_check))

iDEA_to_plot2 <- iDEA_to_plot2 %>%
  mutate(condition = as.character(condition)) %>% 
  mutate(type_of_samples = as.character(type_of_samples)) %>% 
  mutate(Gene.names = as.character(Gene.names))  

iDEA_to_plot2 <- iDEA_to_plot2 %>%
  mutate(condition = factor(condition, levels = c("Downregulated","Non-stablepair","Non-Dysregulated","Upregulated"))) %>%
  mutate(type_of_samples = factor(type_of_samples, levels = c("Cohort","Problem Control","Problem Illness"))) %>% 
  mutate(Gene.names = factor(Gene.names, levels = c(iDEA_genes))) %>% 
  arrange(condition,Gene.names,type_of_samples)

iDEA_to_plot2$combined <- factor(paste(iDEA_to_plot2$condition, 
                                       iDEA_to_plot2$Gene.names, 
                                       iDEA_to_plot2$type_of_samples, sep = "_"), 
                                 levels = unique(paste(iDEA_to_plot2$condition, 
                                                       iDEA_to_plot2$Gene.names, 
                                                       iDEA_to_plot2$type_of_samples, sep = "_")))
iDEA_to_plot2 <- iDEA_to_plot2 %>% 
  mutate(Stripe = ifelse(condition == "Non-stablepair","stripe","none"))

## Barplot
REOA_bars2 <- ggplot(iDEA_to_plot2, aes(y=MED, x=combined))+
  scale_x_discrete(limits = levels(iDEA_to_plot2$combined),
                   labels = iDEA_to_plot2$Gene.names_1)+
  geom_bar(mapping = aes(fill = type_of_samples, color = condition), 
           position="dodge", stat="identity", linewidth = 1)+
  scale_fill_manual(values = type_of_samples)+
  scale_color_manual(values = dynamics_)+
  # geom_bar_pattern(fill = "white",alpha = 0, aes(pattern = Stripe, color = condition),
  #                  position = "dodge", stat = "identity",linewidth=1)+
  # scale_pattern_identity()+
  theme_minimal() +
  labs(title = "REOA results\nTop5 Up and Down proteins",
       x = "Protein")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(size = "none")

REOA_bars2
ggsave(REOA_bars2, filename = "./plots/3_med_sample/REOA_barmap_top_dysregulated.pdf")


#### Pathway analysis on up and downregulated Proteins
## Subset the up and downregulated genes
up_prots <- unique(iDEA_to_plot$Gene.names_1[iDEA_to_plot$condition == "Upregulated"])
down_prots <- unique(iDEA_to_plot$Gene.names_1[iDEA_to_plot$condition == "Downregulated"])

## Define the backGenes
# Here I define the backGenes as the detected genes in a least a sample
backGenes <- unique(iDEA_to_plot$Gene.names_1[iDEA_to_plot$condition != "Non-stablepair"])

## Do the analysis
# Up 
BP_up <- clusterProfiler::enrichGO(gene = up_prots, universe = backGenes,
                                     OrgDb = "org.Hs.eg.db",ont = "BP",
                                     pvalueCutoff = 0.05,qvalueCutoff = 0.1, minGSSize = 5, 
                                     keyType = "SYMBOL")
CC_up <- clusterProfiler::enrichGO(gene = up_prots, universe = backGenes,
                                   OrgDb = "org.Hs.eg.db",ont = "CC",
                                   pvalueCutoff = 0.05,qvalueCutoff = 0.1, minGSSize = 5, 
                                   keyType = "SYMBOL")
MF_up <- clusterProfiler::enrichGO(gene = up_prots, universe = backGenes,
                                   OrgDb = "org.Hs.eg.db",ont = "MF",
                                   pvalueCutoff = 0.05,qvalueCutoff = 0.1, minGSSize = 5, 
                                   keyType = "SYMBOL")
# Down
BP_down <- clusterProfiler::enrichGO(gene = down_prots, universe = backGenes,
                                   OrgDb = "org.Hs.eg.db",ont = "BP",
                                   pvalueCutoff = 0.05,qvalueCutoff = 0.1, minGSSize = 5, 
                                   keyType = "SYMBOL")
CC_down <- clusterProfiler::enrichGO(gene = down_prots, universe = backGenes,
                                   OrgDb = "org.Hs.eg.db",ont = "CC",
                                   pvalueCutoff = 0.05,qvalueCutoff = 0.1, minGSSize = 5, 
                                   keyType = "SYMBOL")
MF_down <- clusterProfiler::enrichGO(gene = down_prots, universe = backGenes,
                                   OrgDb = "org.Hs.eg.db",ont = "MF",
                                   pvalueCutoff = 0.05,qvalueCutoff = 0.1, minGSSize = 5, 
                                   keyType = "SYMBOL")

## Filter the tables by minsize
BP_up <- clusterProfiler::gsfilter(BP_up, by = "Count", min = 5)
CC_up <- clusterProfiler::gsfilter(CC_up, by = "Count", min = 5)
MF_up <- clusterProfiler::gsfilter(MF_up, by = "Count", min = 5)

BP_down <- clusterProfiler::gsfilter(BP_down, by = "Count", min = 5)
CC_down <- clusterProfiler::gsfilter(CC_down, by = "Count", min = 5)
MF_down <- clusterProfiler::gsfilter(MF_down, by = "Count", min = 5)


## Result tables 
# Up 
BP_up_table <- BP_up@result %>% 
  filter(p.adjust <= 0.05)
CC_up_table <- CC_up@result %>% 
  filter(p.adjust <= 0.05)
MF_up_table <- MF_up@result %>% 
  filter(p.adjust <= 0.05)

# Down
BP_down_table <- BP_down@result %>% 
  filter(p.adjust <= 0.05)
CC_down_table <- CC_down@result %>% 
  filter(p.adjust <= 0.05)
MF_down_table <- MF_down@result %>% 
  filter(p.adjust <= 0.05)


## Generate the elements for the loops
# lists
up_list <- list("BP_up" = BP_up,
                "CC_up" = CC_up,
                "MF_up" = MF_up)

down_list <- list("BP_down" = BP_down,
                  "CC_down" = CC_down,
                  "MF_down" = MF_down)
# genes
up_ges <- up_prots
down_ges <- down_prots

# ontologies
ontologies <- c("BP","CC","MF")

### loops for collapsing
# upreg genes collaspe
for(k in 1:length(up_list)){
  # initial elements of the loop
  res <- up_list[[k]]
  ont <- ontologies[k]
  genes <- up_ges
  # collapse with the loop
  if (nrow(res) == 0){
    next
  } else {
  collapes_res <- collapseGO(functional_annot = res@result,
                             pathways = res@geneSets,
                             genes = genes, mingsize = 5, 
                             ontology_to_look = ontologies[k])
  # extract the main functions and save the result as file
  main_functions <- res@result[res@result$ID %in% collapes_res$mainPaths,]
  openxlsx::write.xlsx(x = main_functions, file = paste0("./results/functional_annotation/3_med_sample/",names(up_list)[k],".xlsx"))
  # # do the plot
  go1_collapsed <- ggplot(top_n(main_functions,
                                n = nrow(main_functions), wt = -p.adjust),
                          aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
    geom_bar(stat = "identity", width = 0.5)+
    ggtitle(paste0(names(up_list)[k]))+
    ylab("")+xlab("")+theme_bw()+ labs(fill='FDR') +
    scale_fill_continuous(labels = scientific_format()) ##+
    # theme(legend.text = element_text(size = 20, face = "bold"),
    #       axis.ticks = element_line(colour = "black", size = 1),
    #       panel.border = element_rect(colour = "black", size = 1, fill = NA),
    #       axis.text = element_text(size = 15, face = "bold"),
    #       legend.key.height = unit(1.8,"cm"),
    #       legend.key.width = unit(1.8,"cm"),
    #       legend.title = element_blank())

  ggsave(plot  =go1_collapsed, filename = paste0("./plots/3_med_sample/",names(up_list)[k],".pdf"),
        width = 16, height = 4.5, units = "cm")}
  # pdf(paste0("./plots/",names(up_list)[k],".pdf"), width = 16, height = 4.5, bg = NULL
  res_original <- res@result %>% 
    filter(p.adjust <= 0.05)
  soup <- minestrone(collapsed_list = collapes_res, original_data = res_original)
  res_original <- res_original %>% 
    mutate(collapsed_to = soup$parent[match(Description,
                                            soup$child)])
  res_original <- res_original %>% 
    arrange(collapsed_to,p.adjust)
  openxlsx::write.xlsx(x = res_original, 
                       file = paste0("./results/functional_annotation/3_med_sample/raw_",names(up_list)[k],".xlsx"))
}
# downreg genes collapse
for(k in 1:length(down_list)){
  # initial elements of the loop
  res <- down_list[[k]]
  ont <- ontologies[k]
  genes <- down_ges
  # collapse with the loop
  if (nrow(res) == 0){
    next
  } else {
    collapes_res <- collapseGO(functional_annot = res@result,
                               pathways = res@geneSets,
                               genes = genes, mingsize = 5, 
                               ontology_to_look = ontologies[k])
    # extract the main functions and save the result as file
    main_functions <- res@result[res@result$ID %in% collapes_res$mainPaths,]
    openxlsx::write.xlsx(x = main_functions, file = paste0("./results/functional_annotation/3_med_sample/",names(down_list)[k],".xlsx"))
    # # do the plot
    go1_collapsed <- ggplot(top_n(main_functions,
                                  n = nrow(main_functions), wt = -p.adjust),
                            aes(x=Count, y=reorder(Description,Count), fill=p.adjust))+
      geom_bar(stat = "identity", width = 0.5)+
      ggtitle(paste0(names(down_list)[k]))+
      ylab("")+xlab("")+theme_bw()+ labs(fill='FDR') +
      scale_fill_continuous(labels = scientific_format()) ##+
    # theme(legend.text = element_text(size = 20, face = "bold"),
    #       axis.ticks = element_line(colour = "black", size = 1),
    #       panel.border = element_rect(colour = "black", size = 1, fill = NA),
    #       axis.text = element_text(size = 15, face = "bold"),
    #       legend.key.height = unit(1.8,"cm"),
    #       legend.key.width = unit(1.8,"cm"),
    #       legend.title = element_blank())
    
    ggsave(plot  =go1_collapsed, filename = paste0("./plots/3_med_sample/",names(down_list)[k],".pdf"),
           width = 16, height = 4.5, units = "cm")}
  # pdf(paste0("./plots/",names(down_list)[k],".pdf"), width = 16, height = 4.5, bg = NULL )
  res_original <- res@result %>% 
    filter(p.adjust <= 0.05)
  soup <- minestrone(collapsed_list = collapes_res, original_data = res_original)
  res_original <- res_original %>% 
    mutate(collapsed_to = soup$parent[match(Description,
                                            soup$child)])
  res_original <- res_original %>% 
    arrange(collapsed_to,p.adjust)
  openxlsx::write.xlsx(x = res_original, 
                       file = paste0("./results/functional_annotation/3_med_sample/raw_",names(down_list)[k],".xlsx"))
  
}

### Barplot on the Complement genes
## Get the data
iDEA_to_plot3 <- iDEA_to_plot
selected_rows <- grepl("Complement|complement",iDEA_to_plot3$Protein.Description)
iDEA_to_plot3 <- iDEA_to_plot3[selected_rows,]

iDEA_to_plot3 <- iDEA_to_plot3 %>%
  mutate(condition = as.character(condition)) %>%
  mutate(type_of_samples = as.character(type_of_samples)) %>%
  mutate(Gene.names = as.character(Gene.names))

iDEA_to_plot3 <- iDEA_to_plot3 %>%
  mutate(condition = factor(condition, levels = c("Downregulated","Non-stablepair","Non-Dysregulated","Upregulated"))) %>%
  mutate(type_of_samples = factor(type_of_samples, levels = c("Cohort","Problem Control","Problem Illness"))) %>%
  mutate(Gene.names = factor(Gene.names, levels = c(iDEA_genes))) %>%
  arrange(condition,Gene.names,type_of_samples)

iDEA_to_plot3$combined <- factor(paste(iDEA_to_plot3$condition,
                                       iDEA_to_plot3$Gene.names,
                                       iDEA_to_plot3$type_of_samples, sep = "_"),
                                 levels = unique(paste(iDEA_to_plot3$condition,
                                                       iDEA_to_plot3$Gene.names,
                                                       iDEA_to_plot3$type_of_samples, sep = "_")))
iDEA_to_plot3 <- iDEA_to_plot3 %>%
  mutate(Stripe = ifelse(condition == "Non-stablepair","stripe","none"))
# 
# ## Barplot
REOA_bars3 <- ggplot(iDEA_to_plot3, aes(y=MED, x=combined))+
  scale_x_discrete(limits = levels(iDEA_to_plot3$combined),
                   labels = iDEA_to_plot3$Gene.names_1)+
  geom_bar(mapping = aes(fill = type_of_samples, color = condition),
           position="dodge", stat="identity", linewidth = 1)+
  scale_fill_manual(values = type_of_samples)+
  scale_color_manual(values = dynamics_)+
  # geom_bar_pattern(fill = "white",alpha = 0, aes(pattern = Stripe, color = condition),
  #                  position = "dodge", stat = "identity",linewidth=1)+
  # scale_pattern_identity()+
  theme_minimal() +
  labs(title = "REOA results\nTop5 Up and Down proteins",
       x = "Protein")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(size = "none")

REOA_bars3
ggsave(REOA_bars3, filename = "./plots/3_med_sample/REOA_barmap_COMPLEMENT.pdf")


#### Volcano plot
## Get the data
iDEA_to_plot_volcano <- iDEA_to_plot

# Remove 'Problem Control' samples, Non-stable pairs and 
iDEA_to_plot_volcano <- iDEA_to_plot_volcano %>% 
  mutate(condition = as.character(condition)) %>% 
  mutate(type_of_samples = as.character(type_of_samples)) %>% 
  mutate(Gene.names = as.character(Gene.names)) %>% 
  filter(type_of_samples != "Problem Control") %>% 
  filter(!condition %in% c("Non-stablepair","Non-Dysregulated")) 

# Get the Fold change between the cohort and the 'Problem Ill' protein
wide_idea_to_plot_volcano <- iDEA_to_plot_volcano %>% 
  tidyr::spread(type_of_samples, MED)

wide_idea_to_plot_volcano <- wide_idea_to_plot_volcano %>% 
  subset(select = c(Gene.names,`Cohort`,`Problem Illness`)) %>% 
  group_by(Gene.names) %>% 
  mutate(`Cohort` = min(`Cohort`, na.rm = T)) %>% 
  mutate(`Problem Illness` = min(`Problem Illness`, na.rm = T)) %>% 
  ungroup() %>% 
  distinct() %>% 
  mutate(logFC = `Problem Illness` - `Cohort`)

iDEA_to_plot_volcano <- iDEA_to_plot_volcano %>% 
  mutate(Ill_vs_Cohort = wide_idea_to_plot_volcano$logFC[match(Gene.names,
                                                               wide_idea_to_plot_volcano$Gene.names)]) %>% 
  subset(select = c(Protein.Groups,Accession_1,Gene.names,Gene.names_1,REOA_pval,Ill_vs_Cohort)) %>% 
  distinct()

# The 0 PValues
pval_no_0 <- sort(unique(iDEA_to_plot_volcano$REOA_pval))[2] 

iDEA_to_plot_volcano <- iDEA_to_plot_volcano %>% 
  mutate(REOA_pval = ifelse(REOA_pval == 0,pval_no_0,
                            REOA_pval))
iDEA_to_plot_volcano <- iDEA_to_plot_volcano %>% 
  filter(!is.na(Gene.names_1))

## Volcano plot
# data manipulation
iDEA_to_plot_volcano$fold_change <- ifelse(test = iDEA_to_plot_volcano$Ill_vs_Cohort > 0,
                                    yes = 2^iDEA_to_plot_volcano$Ill_vs_Cohort,
                                    no = (-1/2^iDEA_to_plot_volcano$Ill_vs_Cohort))

iDEA_to_plot_volcano <- iDEA_to_plot_volcano %>%
  mutate(gene_type = ifelse(fold_change >= 0.3785116 & REOA_pval <= 0.05, "Upregulated",
                            ifelse(fold_change <= -0.3785116 & REOA_pval <= 0.05,"Downregulated",
                                   "non-enriched")))

# Define some variables
cols <- c("Downregulated" = "red", 
          "Upregulated" = "blue", 
          "non-enriched" = "grey")

# plotting
vl_plot <- ggplot(data = iDEA_to_plot_volcano,
                  aes(x = Ill_vs_Cohort, y = -log10(REOA_pval), label = Gene.names))+ 
  geom_point(aes(x = Ill_vs_Cohort, y = -log10(REOA_pval),
                 colour = (gene_type)), size = 4)+
  scale_colour_manual(values = cols)+
  scale_size_manual(values = c(6,10))+
  geom_text_repel(data = . %>%
                    filter(Ill_vs_Cohort > 0),
                  nudge_x = 0.50,
                  nudge_y = .1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 2)+
  
  geom_text_repel(data = . %>%
                    filter(Ill_vs_Cohort < 0),
                  nudge_x = -0.50,
                  nudge_y = .1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 2)+
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
  xlab("Change in expression")+
  ylab("log10(Adj.Pvalue)")
vl_plot
ggsave(vl_plot, filename = "./plots/3_med_sample/REOA_vlplot.pdf")

#### Heatmap
### Import the data
raw_data_hm <- merge(expression_mat, annotation, 
                     by.x = "row.names", by.y = "Protein.Groups")

### Create the numeric matrix
raw_data_hm  <- raw_data_hm %>% 
  subset(select = -c(Accession_1,Gene.names_1,Protein.Name,Protein.Description,Row.names)) %>% 
  filter(!duplicated(Gene.names)) %>% 
  filter(!is.na(Gene.names))
rownames(raw_data_hm) <- raw_data_hm$Gene.names
raw_data_hm <- raw_data_hm %>% 
  subset(select = -c(Gene.names))
numeric_values <- as.matrix(raw_data_hm)

### zscore the matrix
# numeric_values <- scale(t(numeric_values))
numeric_values <- t(numeric_values)
### Do the first heatmap
ssGSEA_standard_heatmap <- Heatmap(numeric_values,                            ## The matrix to plot
                                   show_column_names = F, show_row_names = T, ## Show rownames and do not show colnames 
                                   row_names_gp = gpar(fontsize = 10),        ## Rownames graphic options
                                   column_names_gp = gpar(fontsize = 6),      ## Colnames graphic options
                                   cluster_rows = T, cluster_columns = T,     ## Do not cluster rows and colums, the order of the matrix will prevail
                                   # right_annotation = row_annot,              ## The annotation of the 
                                   ####column_km = 4, row_km = 5,            #### BORRAR AIXO !!
                                   name = "Zscored protein \n expression",    ## This is the name that will appear in the legend
                                   column_title = "Integrated protein expression \n profile",) ## Column title acts like if it was a general title
ssGSEA_standard_heatmap <- draw(ssGSEA_standard_heatmap, padding = unit(c(25, ## Down padding
                                                                          5, ## Right padding
                                                                          12.5,## Top padding
                                                                          5), ## Left padding
                                                                        "mm")) ## Padding units
### Group row-wise
## patients_grouping dataframe
patients_grouping <- as.data.frame(tibble::tibble(
  patient = c("A10","A102","A103","A105","A107","A108","A112","A113","A114","A115",
              "A116","A13","A25","A56","A67","A73","A83","A94","A98","A99",
              "C100","C101","C117","C2","C20" ,"C30","C32","C33","C37","C41","C45","C48","C54","C57","C58","C60","C62","C63",
              "C64" ,"C65","C71","E109","E111","E118","E29","E39","E55","E61",
              "E68","E72","E76","E82","E89","E90","E95","C110","mean_sample"),
  grouping = c(rep("Cohort",times = 55),
               "C110","mean_sample")))
rownames(patients_grouping) <- patients_grouping$patient

## Grouping the patients
expression_matrix <- numeric_values

result_matrix <- sapply(unique(patients_grouping$grouping), function(group) {
  groups_patients <- patients_grouping$patient[patients_grouping$grouping == group]
  group_expression <- expression_matrix[grepl(paste0("^", groups_patients, "$", collapse = "|"), rownames(expression_matrix)), 
                                        , drop = FALSE]
  group_mean <- colMeans(group_expression, na.rm = TRUE)
  group_mean[is.na(group_mean)] <- group_expression[which(is.na(group_mean), arr.ind = TRUE)]
  return(group_mean)
})

numeric_values <- result_matrix
remove(result_matrix)
numeric_values <- as.data.frame(numeric_values)
numeric_values <- as.matrix(numeric_values)
numeric_values <- t(numeric_values)

## Proteins to plot
prots_to_plot <- unique(iDEA_to_plot_volcano$Gene.names)
numeric_values <- numeric_values[,prots_to_plot]

### Create the annotations  
## Row annotations
group_annot <- as.data.frame(tibble::tibble(
  pats = rownames(numeric_values),
  Groups = rownames(numeric_values),
))
rownames(group_annot) <- group_annot$pats
group_annot <- group_annot[,2, drop = F]
colours_row <- list("Groups" = c("Cohort" = "grey",
                             "C110" = "yellow",
                             "mean_sample" = "purple"))
row_annot <- HeatmapAnnotation(df = group_annot,
                               col = colours_row,
                               which = "row")
## Column annotations
prots_annot <- iDEA_to_plot_volcano
prots_annot <- prots_annot %>% 
  subset(select = c(Gene.names, gene_type)) %>% 
  arrange(desc(gene_type)) %>% 
  distinct()
prots_annot <- as.data.frame(prots_annot)
rownames(prots_annot) <- prots_annot$Gene.names
prots_annot <- prots_annot[,-1, drop = F]
colnames(prots_annot)[1] <- "Direction of expression"
colours_col <- list("Direction of expression" = c("Upregulated" = "red",
                                                  "Downregulated" = "blue"))
col_annot <- HeatmapAnnotation(df = prots_annot,
                               col = colours_col,
                               which = "column")

# Reorder the numeric_values matrix
numeric_values <- numeric_values[rownames(group_annot),
                                 rownames(prots_annot)]

### Create the splits
split_row <- factor(group_annot$Groups, levels =c("C110","mean_sample","Cohort"))
split_col <- factor(prots_annot$`Direction of expression`, levels =c("Upregulated","Downregulated"))

numeric_values <- scale(numeric_values)

###  Do the second heatmap
hmap_rowise <- Heatmap(numeric_values, 
                      cluster_row_slices = FALSE, 
                      row_split = split_row,
                      column_split = split_col,
                      show_row_dend = F, show_column_dend = F, 
                      show_heatmap_legend = T,
                      column_names_gp = gpar(fontsize = 10), 
                      column_title_gp = gpar(fontsize = 12),
                      right_annotation = row_annot, 
                      top_annotation = col_annot,
                      column_title = "REOA diff expression results", 
                      name = "log2 scaled\nintensity",
                      cluster_columns = F, show_row_names = T)

pdf("./plots/3_med_sample/REOA_hmap.pdf", width = 15, height = 10)
hmap_rowise_drawn <- draw(hmap_rowise)
dev.off()





