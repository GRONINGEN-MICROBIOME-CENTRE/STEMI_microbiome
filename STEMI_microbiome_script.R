### Analyses for cardiovascular risk spectrum microbiome
### March 2023 - Femke Prins

#=====Preparing data=====# ----
#### Load libraries ----
library(readr)
library(tidyverse)
library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ape)
library(stringr)
library(ggpubr)
library(reshape)
library(xlsx)
library(plotly)
library(rmarkdown)
library(Maaslin2)
library(Hmisc)
library(corrplot)
library(compositions)
library(table1)
library(foreign)
library(reshape2)
library(grid)
library(patchwork)
library(forcats)
library(ggsignif)
library(emmeans)
library(ggvenn)
#### Importing clean taxa data (participants selected from total cohort, taxa's with only zeros removed) ----
Taxa_STEMI <- read.table(file = "~/Documents/MDPhD/STEMI_microbiome/Taxa_STEMI_clean", header = TRUE, sep = "\t",row.names = 1)
Taxa_DAG3 <- read.table(file = "~/Documents/MDPhD/STEMI_microbiome/Taxa_DAG3_clean", header = TRUE, sep = "\t",row.names = 1)

#### Rescaling it with Ranko's function and selecting bacteria and archaea (no further filtering) ----
filterMetaGenomeDF <- function(inDF,presPerc = 0.1,minMRelAb = 0.01,minMedRelAb=0.0,
                               rescaleTaxa=T,verbose=T,
                               keepDomains=c('Bacteria','Archaea'),
                               keepLevels=c('S','G','F','O','C','P'),
                               keepUnknown=F) {
  
  # -- drop unknown & rescale? --
  if (!keepUnknown) {
    inDF$UNKNOWN <- NULL
  }
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  tColsNMG = grep('k__',colnames(inDF),invert = T) # colums with microbiome
  # replaces NAs in microbiome with 0s
  for (c in tCols) {
    inDF[,c][is.na(inDF[,c])] <- 0.0
  }
  
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    nrnZ = as.numeric(sum(inDF[,c]!=0.0))
    if ( (nrnZ/as.numeric(nrow(inDF))) < presPerc) {
      # if (verbose) {
      #   print (paste('col',c,': ',colnames(inDF)[c],'; nr non Zero:',nrnZ,'=',nrnZ/as.numeric(nrow(inDF)),'>> Column removed!'))
      # }
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('^[dgtspcfko]__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = mean(inDF[,c])
    if ( mn < minMRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; mean rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }
  
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = median(inDF[,c])
    if ( mn < minMedRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; median rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }  
  
  # keep domains
  # -----------------------------
  toKeep <- NULL
  if (length(keepDomains) == 1) {
    if (keepDomains == 'All' | keepDomains == "") {
      toKeep = grep('k__',colnames(inDF),invert = T)
    } else {
      toKeep = grep(paste('k__',d,sep=''),colnames(inDF))
    }
  } else {
    for (d in keepDomains) {
      #d print(d)
      toKeep = c(toKeep,grep(paste('k__',d,sep=''),colnames(inDF)))
    }
  }
  inDF <- inDF[,toKeep]
  
  # remove taxonomic levels
  # -----------------------------
  inDFnonTaxa <- as.data.frame(inDF[,grep('k__',colnames(inDF),invert=T)])
  colnames(inDFnonTaxa) <- colnames(inDF)[grep('k__',colnames(inDF),invert=T)]
  
  inDF2 <- as.data.frame(inDF[,grep('k__',colnames(inDF),invert=F)])
  # pick strains (T)
  taxaTCols <- grep('t__',colnames(inDF2))
  taxaT <- as.data.frame(inDF2[,taxaTCols])
  colnames(taxaT) <- colnames(inDF2)[grep('t__',colnames(inDF2))]
  if (length(taxaTCols) > 0) {inDF2 <- inDF2[,-taxaTCols]}
  # pick species (S)
  taxaSCols <- grep('s__',colnames(inDF2))
  if (length(taxaSCols) > 0) {
    taxaS <- as.data.frame(inDF2[,taxaSCols])
    colnames(taxaS) <- colnames(inDF2)[grep('s__',colnames(inDF2))]
  }
  if (length(taxaSCols) > 0) {inDF2 <- inDF2[,-taxaSCols]}
  # pick genera (G)
  taxaGCols <- grep('g__',colnames(inDF2))
  if (length(taxaGCols) > 0) {
    taxaG <- as.data.frame(inDF2[,taxaGCols])
    colnames(taxaG) <- colnames(inDF2)[grep('g__',colnames(inDF2))]
  }
  if (length(taxaGCols) > 0) {inDF2 <- inDF2[,-taxaGCols]}
  # pick families (F)
  taxaFCols <- grep('f__',colnames(inDF2))
  if (length(taxaFCols) > 0) {    
    taxaF <- as.data.frame(inDF2[,taxaFCols])
    colnames(taxaF) <- colnames(inDF2)[grep('f__',colnames(inDF2))]
  }
  if (length(taxaFCols) > 0) {inDF2 <- inDF2[,-taxaFCols]}
  # pick orders (O)
  taxaOCols <- grep('o__',colnames(inDF2))
  if (length(taxaOCols) > 0) {
    taxaO <- as.data.frame(inDF2[,taxaOCols])
    colnames(taxaO) <- colnames(inDF2)[grep('o__',colnames(inDF2))]
  }
  if (length(taxaOCols) > 0) {inDF2 <- inDF2[,-taxaOCols]}
  # pick classes (C)
  taxaCCols <- grep('c__',colnames(inDF2))
  if (length(taxaCCols) > 0) {
    taxaC <- as.data.frame(inDF2[,taxaCCols])
    colnames(taxaC) <- colnames(inDF2)[grep('c__',colnames(inDF2))]
  }
  if (length(taxaCCols) > 0) {
    if( length(colnames(inDF2)) - length(taxaCCols) == 1) {
      ccn <- colnames(inDF2)[grep('c__',colnames(inDF2))]
      tempN <- colnames(inDF2)[!(colnames(inDF2) %in% ccn)]
      inDF2 <- as.data.frame(inDF2[,-taxaCCols])
      colnames(inDF2) <- tempN
    } else {
      inDF2 <- inDF2[,-taxaCCols]
    }
  }
  # pick phyla (P)
  taxaPColsKeepN <- NULL
  taxaPCols <- grep('p__',colnames(inDF2))
  if (length(taxaPCols) > 0) {
    taxaPColsN <- colnames(inDF2)[grep('p__',colnames(inDF2))]
    taxaPColsKeep <- grep('p__',colnames(inDF2),invert = T)
    taxaPColsKeepN <- colnames(inDF2)[grep('p__',colnames(inDF2),invert = T)]
    taxaP <- as.data.frame(inDF2[,taxaPCols])
    colnames(taxaP) <- taxaPColsN
  }
  # pick kingdoms (K)
  if (length(taxaPCols) > 0) {inDF2 <- as.data.frame(inDF2[,taxaPColsKeep])}
  if (!is.null(taxaPColsKeepN)){
    colnames(inDF2) <- taxaPColsKeepN
  }
  taxaK <- inDF2
  # pick 
  oDF <- inDFnonTaxa
  if (verbose) {print ('Keeping following taxonomic levels:'); print(keepLevels)}
  if (ncol(taxaK) > 0) {  
    if ('K' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaK),'Kingdoms'))}
      if (rescaleTaxa) {taxaK <- taxaK/rowSums(taxaK)}
      oDF <- cbind(oDF,taxaK)}
  }
  if (ncol(taxaP) > 0) {
    if ('P' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaP),'Phyla'))} 
      if (rescaleTaxa) {taxaP <- taxaP/rowSums(taxaP)}
      oDF <- cbind(oDF,taxaP)}
  }
  if (ncol(taxaC) > 0) {
    if ('C' %in% keepLevels) {
      if (verbose) {print(paste(' -> kept',ncol(taxaC),'Classes'))} 
      if (rescaleTaxa) {taxaC <- taxaC/rowSums(taxaC)}
      oDF <- cbind(oDF,taxaC)}
  }
  if (ncol(taxaO) > 0) {
    if ('O' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaO),'Orders'))}
      if (rescaleTaxa) {taxaO <- taxaO/rowSums(taxaO)}
      oDF <- cbind(oDF,taxaO)}
  }
  if (ncol(taxaF) > 0) {
    if ('F' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaF),'Families'))}
      if (rescaleTaxa) {taxaF <- taxaF/rowSums(taxaF)}
      oDF <- cbind(oDF,taxaF)}
  }
  if (ncol(taxaG) > 0) { 
    if ('G' %in% keepLevels) {if (verbose){ print(paste(' -> kept',ncol(taxaG),'Genera'))}
      if (rescaleTaxa) {taxaG <- taxaG/rowSums(taxaG)}
      oDF <- cbind(oDF,taxaG)}
  }
  if (ncol(taxaS) > 0) {
    if ('S' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaS),'Species'))}
      if (rescaleTaxa) {taxaS <- taxaS/rowSums(taxaS)}
      oDF <- cbind(oDF,taxaS)}
  }
  if (ncol(taxaT) > 0) {
    if ('T' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaT),'Strains'))}
      if (rescaleTaxa) {taxaT <- taxaT/rowSums(taxaT)}
      oDF <- cbind(oDF,taxaT)}
  }
  if (verbose) {print ('data processing done, returning Dataframe')}
  oDF
}

Taxa_DAG3_rel <- filterMetaGenomeDF(Taxa_DAG3,presPerc = 0.0,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                    keepDomains=c('Bacteria','Archaea'),keepLevels=c('S','G','F','O','C','P','K'),
                                    keepUnknown=F) #keeps 676 variables

Taxa_STEMI_rel <- filterMetaGenomeDF(Taxa_STEMI,presPerc = 0.0,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                     keepDomains=c('Bacteria','Archaea'),keepLevels=c('S','G','F','O','C','P','K'),
                                     keepUnknown=F) #keeps 475 variables
#### Merging datasets ----
#Preparing for merging (microbiome)
Taxa_DAG3_rel <- rownames_to_column(Taxa_DAG3_rel, var = "ID")
Taxa_STEMI_rel <- rownames_to_column(Taxa_STEMI_rel, var = "ID")
Taxa_merged_all <- full_join(Taxa_DAG3_rel, Taxa_STEMI_rel) #keeps 706 variables
rownames(Taxa_merged_all) = Taxa_merged_all$ID # making ID the rowname again
Taxa_merged_all[,1] <- NULL

#Removing taxa that are not in both cohorts
Taxa_merged_all %>% select_if(~ !any(is.na(.))) -> Taxa_merged #keeps 445 variables

#For species, keep s__Ruminococcus_sp_CAG_488 because of prevalence >5% in cohort DAG3
Taxa_merged_all %>% select(c("k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruminococcus.s__Ruminococcus_sp_CAG_488")) -> Taxa_merge_Ruminococcus
Taxa_merged <- merge(Taxa_merged, Taxa_merge_Ruminococcus, by = 0)
rownames(Taxa_merged) = Taxa_merged$Row.names
Taxa_merged[,1] <- NULL
Taxa_merged %>% mutate_at(c(446), ~replace(., is.na(.), 0)) -> Taxa_merged #replace NA in STEMI cohort to zero's


#### Importing metadata (see script Merge_metadata.R) ----
Meta_merged<- readRDS(file = "~/Documents/MDPhD/STEMI_microbiome/Merged_metadata.Rds")
rownames(Meta_merged) <- Meta_merged$sample_ID
Meta_merged %>% dplyr::select(c(1:12, 14:35, 49:50)) -> Meta_merged #keep relevant metadata

#### Grabbing taxonomic level ('S','G','F','O','C','P') ----
Taxa_merged_S <- Taxa_merged[,grepl('s__',colnames(Taxa_merged),ignore.case = T)] #keeps 240 species

Taxa_merged_G <- Taxa_merged[,grepl('g__',colnames(Taxa_merged),ignore.case = T)] #keeps 101 genera
Taxa_merged_G<- Taxa_merged_G[,!grepl('s__',colnames(Taxa_merged_G),ignore.case = T)]

Taxa_merged_F <- Taxa_merged[,grepl('f__',colnames(Taxa_merged),ignore.case = T)] #keeps 45 families
Taxa_merged_F<- Taxa_merged_F[,!grepl('g__',colnames(Taxa_merged_F),ignore.case = T)]

Taxa_merged_O <- Taxa_merged[,grepl('o__',colnames(Taxa_merged),ignore.case = T)] #keeps 28 orders
Taxa_merged_O<- Taxa_merged_O[,!grepl('f__',colnames(Taxa_merged_O),ignore.case = T)]

Taxa_merged_C <- Taxa_merged[,grepl('c__',colnames(Taxa_merged),ignore.case = T)] #keeps 20 classes
Taxa_merged_C<- Taxa_merged_C[,!grepl('o__',colnames(Taxa_merged_C),ignore.case = T)]

Taxa_merged_P <- Taxa_merged[,grepl('p__',colnames(Taxa_merged),ignore.case = T)] #keeps 10 phyla
Taxa_merged_P<- Taxa_merged_P[,!grepl('c__',colnames(Taxa_merged_P),ignore.case = T)]


#=====Analyses:species=====#----
#### Alpha diversity ----
calculate_alpha <- function(inDF,IDcol="RowNames",
                            metrics=c("shannon","simpson","invsimpson","richness"),
                            DIVlvls=c("taxS")) {
  # select IDs column
  if (IDcol == "RowNames") {
    DIVMatrix <- data.frame(RN=rownames(inDF))
  } else {
    DIVMatrix <- data.frame(IDcol=inDF[[IDcol]])
    colnames(DIVMatrix)[1] <- IDcol
  }
  # iterate over metrics, calculate each
  # NOTE: richness is not implemented in vegan, requires special treatment
  for (l in DIVlvls) {
    if (grepl('^tax.$',l)) {
      toUse <- gsub('^tax','',l) } }
  for (m in metrics) {
    print(paste0('  > calculating ',m,'[',l,']'))
    if (m=="richness") {
      inDFpa <- inDF
      inDFpa[inDFpa > 0] <- 1
      dv <- rowSums(inDFpa)
    } else {
      dv <- diversity(inDF,index = m)
    }
    DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
    colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.',toUse,'.',m)
  }
  return(DIVMatrix)
}
Alpha_metrices <- calculate_alpha(Taxa_merged_S, IDcol="RowNames", metrics=c("shannon","simpson","invsimpson","richness"), DIVlvls=c("taxS"))
#Alpha_metrices <- calculate_alpha(Taxa_CLR_S, IDcol="RowNames", metrics=c("shannon","simpson","invsimpson","richness"), DIVlvls=c("taxS"))
Data_div_total <- merge(Alpha_metrices, Meta_merged, by = 0)
rownames(Data_div_total) <- Data_div_total$Row.names
Data_div_total$Row.names <- NULL
Data_div_total$richness_reads = Data_div_total$DIV.S.richness/Data_div_total$reads

#Plotting violin plots
my_comparisons_all=list(c("STEMI","Low"),c("STEMI","Intermediate"),c("STEMI","High"),c("Low","Intermediate"),c("Low","High"),c("Intermediate","High"))
my_cols <- c("STEMI"="#b91a1a", "Low"="#00938d", "Intermediate"="#36528f", "High"="#4b0257")
Data_div_total$Risk_group <- ordered(Data_div_total$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI"))
Data_div_total$Risk_group <- factor(Data_div_total$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI"))


a_shannon <- Data_div_total %>% 
  ggplot(aes(x=Risk_group, y=DIV.S.shannon, fill=Risk_group)) +
  geom_jitter(aes(color=Risk_group), alpha=0.5, width=0.1)+
  geom_boxplot(aes(color=Risk_group), alpha=0, width=0.3, size=0.8)+
  geom_violin(aes(color = Risk_group), trim = FALSE, position = position_dodge(0.9), alpha=0, size=0.8)+
  theme_classic()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
  theme(legend.position = "none")+
  ylab("Shannon Diversity")+
  scale_color_manual(values = my_cols)+
  stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", hide.ns = FALSE)
a_shannon

a_simpson <- Data_div_total %>% 
  ggplot(aes(x=Risk_group, y=DIV.S.simpson, fill=Risk_group)) +
  geom_jitter(aes(color=Risk_group), alpha=0.5, width=0.1)+
  geom_boxplot(aes(color=Risk_group), alpha=0, width=0.3, size=0.8)+
  geom_violin(aes(color = Risk_group), trim = FALSE, position = position_dodge(0.9), alpha=0, size=0.8)+
  theme_classic()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
  theme(legend.position = "none")+
  ylab("Simpson Diversity")+
  scale_color_manual(values = my_cols)+
  stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", hide.ns = FALSE)
a_simpson

a_invsimpson <- Data_div_total %>% 
  ggplot(aes(x=Risk_group, y=DIV.S.invsimpson, fill=Risk_group)) +
  geom_jitter(aes(color=Risk_group), alpha=0.5, width=0.1)+
  geom_boxplot(aes(color=Risk_group), alpha=0, width=0.3, size=0.8)+
  geom_violin(aes(color = Risk_group), trim = FALSE, position = position_dodge(0.9), alpha=0, size=0.8)+
  theme_classic()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
  theme(legend.position = "none")+
  ylab("Inverse Simpson Diversity")+
  scale_color_manual(values = my_cols)+
  stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", hide.ns = FALSE)
a_simpson

a_richness <- Data_div_total %>% 
  ggplot(aes(x=Risk_group, y=DIV.S.richness, fill=Risk_group)) +
  geom_jitter(aes(color=Risk_group), alpha=0.5, width=0.1)+
  geom_boxplot(aes(color=Risk_group), alpha=0, width=0.3, size=0.8)+
  geom_violin(aes(color = Risk_group), trim = FALSE, position = position_dodge(0.9), alpha=0, size=0.8)+
  theme_classic()+
  theme(
    axis.title.y = element_text(size=14, face="bold", colour = "black"),    
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black", hjust = 0),
    legend.text = element_text(face="bold"),
    legend.title = element_text(face = "bold"),
    plot.title=element_text(face = "bold"))+
  theme(legend.position = "none")+
  ylab("Richness")+
  scale_color_manual(values = my_cols)+
  stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", hide.ns = FALSE)
a_richness

#Plotting all diversities in one plot
Combined_a_diversities <- ggarrange(a_shannon, a_simpson, a_invsimpson, a_richness, ncol = 2, nrow = 2)
print(Combined_a_diversities)

#calculate means
Data_div_total %>% 
  group_by(Risk_group) %>% 
  summarise_at(vars(DIV.S.richness), list(Richness = ~mean(.), SD = ~sd(.)))

#### CLR transformation (after filtering on prevalence) ----
#First filtering on prevalence (10%)
Prevalence.filter <- function(x, threshold = 10){
  remove_cols <- vector()
  for (i in 1:ncol(x)) {
    cname <- colnames(x)[i]
    if(colSums(x[i] != 0) / nrow(x) *100 < threshold){
      remove_cols <- c(remove_cols, cname) 
    }
  }
  x_new <- x %>% dplyr::select(-remove_cols)
  print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))
  return(x_new)
}
taxa_S_filter <- Prevalence.filter(Taxa_merged_S, threshold = 10)

#Do CLR transformation
my_pseudocount=min(taxa_S_filter[taxa_S_filter!=0])/2.
Taxa_CLR_S <- decostand(taxa_S_filter, method = "clr", pseudocount = my_pseudocount)


#### Beta diversity calculations and plotting ----
#Calculating beta diversity
vegdist(Taxa_CLR_S, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by = 'row.names')

#Making plots
PCoA_meta$Risk_group <- ordered(PCoA_meta$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI"))

ref <- reformulate("Risk_group","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids
centroids_long <- gather(centroids, coordinate, mean, V1:V5, factor_key=TRUE)
write.xlsx(centroids_long,"~/Documents/MDPhD/Results_STEMI_DAG3/coordinates_mean.xlsx")

my_cols <- c("STEMI"="#b91a1a", "Low"="#00938d", "Intermediate"="#36528f", "High"="#4b0257")

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=Risk_group)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=Risk_group)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=Risk_group)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=Risk_group)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

#test differences between centroids
wilcoxVarsTouse <- c("V1","V2","V3","V4","V5")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:5]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, Risk_group)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$Risk_group, p.adjust.method="none") -> wilcox_df
  melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

Wilcox_PCoA_Results <- transform(Wilcox_PCoA_Results, combination=paste(Variable1, Variable2, sep=":")) #merge Group1, Group2 variables into a new variable
subset(Wilcox_PCoA_Results, !is.na(P_value)) -> Wilcox_PCoA_Results_2
Wilcox_PCoA_Results_2 %>% select(c(Coordinate, P_value, combination)) -> Wilcox_PCoA_Results_2
Wilcox_PCoA_Results_2_wide <- spread(Wilcox_PCoA_Results_2, combination, P_value)
write.xlsx(Wilcox_PCoA_Results_2_wide,"~/Documents/MDPhD/Results_STEMI_DAG3/coordinates_pvalues.xlsx")

#Betadisper
Merged_all <- merge(Taxa_CLR_S, Meta_merged, by = 'row.names')
rownames(Merged_all) <- Merged_all$Row.names
Merged_all$Row.names <- NULL
colnames(Merged_all)
dis <- vegdist(Merged_all[,1:136],method = "euclidean") #create distance matrix
#calculate multivariate dispersions
mod <- betadisper(dis, Merged_all$Risk_group)
anova(mod)
TukeyHSD(mod) -> mod_HSD
plot(mod_HSD)
permutest(mod, pairwise = TRUE, permutations = 99) -> pmod

pstat <- permustats(pmod)
densityplot(pstat, scales = list(x = list(relation = "free")))
qqmath(pstat, scales = list(relation = "free"))

#plot the findings
cols_colours <- c("#00938d","#36528f","#4b0257","#b91a1a")
mod$group <- ordered(mod$group, levels = c("Low", "Intermediate", "High", "STEMI"))
boxplot(mod, xlab = "", col = cols_colours)

#### Adonis (univariate) ----
# load phenotypes
inPhenos <- Meta_merged
# make character variables factors
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)
ordered(inPhenos$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI")) -> inPhenos$Risk_group
str(inPhenos) #check

# remove IDs and other (irrelevant) phenotypes
adonisVarsTouse <- colnames(inPhenos)
toExclude = c("sample_ID", "sex", "age", "sys_bp", "dys_bp", "diabetes", "oral_antidiabetic", "insulin", "hdl", "cholesterol_total", "ldl",
              "smoking_status", "at_antagonist", "beta_block", "ca_antagonist", "ace_inhibitor","diuretics_low_ceiling",
              "hypercholesterolemia", "hypertension", "triglyceride", "Framingham")
adonisVarsTouse <- adonisVarsTouse[!adonisVarsTouse %in% toExclude]

adonisResults <- NULL
permNR <- 1000
inMB <- Taxa_CLR_S
inMB <- tibble::rownames_to_column(inMB, "sample_ID")

for (i in adonisVarsTouse[1:15]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"sample_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="sample_ID",y=inMB,by.y="sample_ID")
  rownames(allDF) <- allDF$sample_ID
  allDF$sample_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis(inBC ~ av,permutations=permNR)
    aov_table <- ad$aov.tab
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov_table[1,1],
                         SumsOfSqs=aov_table[1,2],
                         MeanSqs=aov_table[1,3],
                         FModel=aov_table[1,4],
                         R2=aov_table[1,5],
                         pval=aov_table[1,6],
                         FDR.BH=NA,
                         Significant=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    write.table(adonisResults,paste0("~/Documents/MDPhD/Results_STEMI_DAG3/adonis_results_table.csv"),sep=",",row.names=F)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$Var <- NULL
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant="No"
adonisResults$Significant[adonisResults$FDR.BH<0.05]="Yes"

write.xlsx(adonisResults,"~/Documents/MDPhD/Results_STEMI_DAG3/adonis_univariate.xlsx")


#univariate alle phenotypes testen, goed voor het verhaal

my_colours <- c("Yes"="darkolivegreen3", "No"="indianred3")
#plot for species
adonis_species <- ggplot(adonisResults, aes(reorder(row.names(adonisResults), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Taxonomy = Species (univariate)") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = my_colours)
plot(adonis_species)

#### Adonis (multivariate) ----
All_df <- merge(Meta_merged, Taxa_CLR_S, by = 'row.names')
ad_meta <- All_df[c(1:37)]
ad_taxa <- All_df[c(38:173)]
ad_taxa_dm <- vegdist(ad_taxa,method = "euclidean")
ad1 <- adonis2(formula = ad_taxa_dm ~ proton_pump_inhib + Risk_group + reads, data = ad_meta, permutations = 1000, by = "margin")

model.frame(ad1) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"

adonis_species_multi <- ggplot(adonis_mv, aes(reorder(row.names(adonis_mv), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Taxonomy = Species (multivariate)") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = my_colours)
plot(adonis_species_multi)

write.xlsx(adonis_mv,"~/Documents/MDPhD/Results_STEMI_DAG3/adonis_multivariate.xlsx")


#Risk group seems to be significant

#### Adonis with risk score variables ----
#Which components or risk score have biggest effect, re-run same model with covariates
All_df <- merge(Meta_merged, Taxa_CLR_S, by = 'row.names')
rownames(All_df) <- All_df$Row.names
All_df %>% mutate(medication_hypertension= ifelse(diuretics_low_ceiling=="Y" | beta_block=="Y" | ca_antagonist=="Y" | ace_inhibitor =="Y" | at_antagonist=="Y", "Y", "N" ) ) -> All_df
ad_meta <- All_df[c(1:37, 174)]
ad_taxa <- All_df[c(38:173)]
ad_taxa_dm <- vegdist(ad_taxa,method = "euclidean")
ad2 <- adonis2(formula = ad_taxa_dm ~ sex + age + sys_bp + diabetes + hdl + cholesterol_total + medication_hypertension + 
                 smoking_status + reads + proton_pump_inhib, data = ad_meta, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv_risk
adonis_mv_risk$Significant="No"
adonis_mv_risk$Significant[adonis_mv_risk$`Pr(>F)`<0.05]="Yes"

adonis_species_multi_riskscore <- ggplot(adonis_mv_risk, aes(reorder(row.names(adonis_mv_risk), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Taxonomy = Species. Riskscore variables (multivariate)") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = my_colours)
plot(adonis_species_multi_riskscore)

write.xlsx(adonis_mv_risk,"~/Documents/MDPhD/Results_STEMI_DAG3/adonis_riskscore_variables.xlsx")

#### Differential abundance (emmeans), FDR correction per contrast ----
df_lm1 <- Meta_merged[c("sample_ID", "proton_pump_inhib", "reads", "Risk_group")]
df_lm1[sapply(df_lm1, is.character)] <- lapply(df_lm1[sapply(df_lm1, is.character)], as.factor)
df_lm1$sample_ID <- as.character(df_lm1$sample_ID)
df_lm1$Risk_group <- factor(df_lm1$Risk_group, ordered = FALSE)
df_lm1$Risk_group <- factor(df_lm1$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI"))
levels(df_lm1$Risk_group)
Taxa_CLR_S$sample_ID <- rownames(Taxa_CLR_S)

df_lm1 <- left_join(df_lm1, Taxa_CLR_S, by="sample_ID")
Taxa_CLR_S$sample_ID <- NULL

myContr1 <-list(Low_Intermediate=c(-1,1,0,0),
                Low_High=c(-1,0,1,0),
                Low_STEMI=c(-1,0,0,1),
                Intermediate_High=c(0,-1,1,0),
                Intermediate_STEMI=c(0,-1,0,1),
                High_STEMI=c(0,0,-1,1),
                LIH_STEMI=c(-(1/3),-(1/3),-(1/3),1))

lm_contrasts <- vector("list", ncol(Taxa_CLR_S))
names(lm_contrasts) <- colnames(Taxa_CLR_S)
lm_models <- vector("list", ncol(Taxa_CLR_S))
names(lm_models) <- colnames(Taxa_CLR_S)

for(i in colnames(df_lm1)[grepl(x=colnames(df_lm1), pattern="k__")]) {
  #tryCatch({
  print(i)
  m_species <- lm(df_lm1[,i] ~ Risk_group + proton_pump_inhib + reads, 
                  data=df_lm1[,c(i,"Risk_group", "proton_pump_inhib", "reads")])
  lm_models[[i]] <- m_species
  mm_species <- emmeans::emmeans(m_species,"Risk_group")
  lm_contrasts[[i]] <- as.data.frame(emmeans::contrast(mm_species, myContr1)) } #transform the emmGrid object into a data frame, which will be recognized as a vector in tibble

#Get all results in one dataframe
lm_raw_results <- lm_contrasts %>% bind_rows(.id="taxon_id") 

#Calculate FDR per contrast
lm_contrasts_sigFDR <- vector("list",7)
names(lm_contrasts_sigFDR) <- c("Low_Intermediate", "Low_High","Low_STEMI","Intermediate_High","Intermediate_STEMI","High_STEMI","LIH_STEMI")
lm_contrasts_allFDR <- vector("list",7)
names(lm_contrasts_allFDR) <- c("Low_Intermediate", "Low_High","Low_STEMI","Intermediate_High","Intermediate_STEMI","High_STEMI","LIH_STEMI")

for(c in c("Low_Intermediate", "Low_High","Low_STEMI","Intermediate_High","Intermediate_STEMI","High_STEMI","LIH_STEMI")) {
  df <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    filter(p_adj < 0.05) #%>%
    #mutate(taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
           # sign=ifelse(estimate>0, "positive", "negative"),
           #contrast=factor(contrast))
  lm_contrasts_sigFDR[[c]] <- df
  
  df1 <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
   mutate(taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
           sign=ifelse(estimate>0, "positive", "negative"),
           contrast=factor(contrast))
  lm_contrasts_allFDR[[c]] <- df1
}

lm_contrasts_sigFDR <- lm_contrasts_sigFDR %>% bind_rows() #significant dataframe
lm_contrasts_allFDR <- lm_contrasts_allFDR %>% bind_rows() #all results
write.xlsx(lm_contrasts_allFDR,"~/Documents/MDPhD/Results_STEMI_DAG3/DDA_species.xlsx")

df_lm1 %>% group_by(Risk_group) %>% summarise_at(vars("k__Bacteria.p__Actinobacteria.c__Coriobacteriia.o__Coriobacteriales.f__Coriobacteriaceae.g__Collinsella.s__Collinsella_stercoris"), list(name = mean))

#### Heatmap plotting ----
complete(lm_contrasts_sigFDR, taxon_id, contrast) -> df_heatmap_species
rng_estimates <- range(df_heatmap_species[!is.na(df_heatmap_species$estimate),]$estimate)

library(stringr)
df_heatmap_species$taxon_id <- str_extract(df_heatmap_species$taxon_id, "s__.+")

df_heatmap_species$taxon_id <- str_replace(df_heatmap_species$taxon_id, "s__", "")
df_heatmap_species$taxon_id <- str_replace(df_heatmap_species$taxon_id, "_", " ")
df_heatmap_species$taxon_id <- str_replace(df_heatmap_species$taxon_id, "_", " ")
#df_heatmap_species$taxon_id <- str_replace(df_heatmap_species$taxon_id, ".Collinsella. massiliensis", "Collinsella massiliensis")
unique(df_heatmap_species$taxon_id)

# Define the order of the levels
desired_order <- c("Low_STEMI", "Intermediate_STEMI", "High_STEMI", "LIH_STEMI", "Intermediate_HS")
# Convert the "contrast" variable to a factor with the desired order of levels
df_heatmap_species$contrast <- factor(df_heatmap_species$contrast, levels = desired_order)

heatmap_species <- df_heatmap_species %>% 
  ggplot(aes(y=taxon_id,x=contrast,fill=estimate)) +
  coord_fixed()+ # make the cell contain square
  geom_tile(colour="black") +
  scale_fill_gradient2(name="Effect size", low = "#1111AA", mid = "white", high = "#AA1111",midpoint = 0, na.value = "grey90",limits=c(floor(rng_estimates[1]), ceiling(rng_estimates[2]))) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.text.y = element_text(size=8, face="italic", colour = "black"),
        axis.text.x = element_text(size=6, face="bold", colour = "black", hjust= 1, angle=45),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.position = "top",
        legend.text = element_text(size=6, face="bold", colour = "black"),
        legend.title =element_text(size=6, face="bold", colour = "black"),
        legend.margin=margin(grid::unit(0, "cm")),
        legend.key.height=grid::unit(0.2, "cm"),
        legend.key.width=grid::unit(0.8, "cm")) +
  labs(x=NULL, y=NULL)
print(heatmap_species)

#### Linear model for trend testing ----
df_lm2 <- df_lm1
df_lm2$Risk_group_num <- as.numeric(df_lm2$Risk_group)
#col_species <- names(df_lm2)
#gsub(".*s__","",col_species) -> col_species
#names(df_lm2) = c(col_species)

lm_all_trend <- data.frame(species = character(), coefficients = numeric(), p_values = numeric(), stringsAsFactors = FALSE)
## all species
for(i in colnames(df_lm2)[5:140]) {
  print(i)
  m_species <- lm(df_lm2[,i] ~ Risk_group_num + proton_pump_inhib + reads, 
                  data=df_lm2[,c(i,"Risk_group_num", "proton_pump_inhib", "reads")])
  # Extract coefficients and p-values
  coef <- coef(m_species)
  p_val <- summary(m_species)$coefficients[, "Pr(>|t|)"]
  # Append species, coefficients, and p-values to the dataframe
  lm_all_trend <- rbind(lm_all_trend, data.frame(species = i, coefficients = coef, p_values = p_val))
}

lm_all_trend$FDR.BH=p.adjust(lm_all_trend$p_values, method = "BH")
lm_all_trend$Significant="No"
lm_all_trend$Significant[lm_all_trend$FDR.BH<0.05]="Yes"
lm_all_trend$coef_name <- rownames(lm_all_trend)
lm_all_trend$p_values <- format(lm_all_trend$p_values, scientific = FALSE)

# Select rows that start with "Risk_group" in coef_name column
selected_rows <- lm_all_trend[grepl("^Risk_group", lm_all_trend$coef_name), ]

write.xlsx(selected_rows,"~/Documents/MDPhD/Results_STEMI_DAG3/linear_trend_species.xlsx")

subset_df <- selected_rows[selected_rows$Significant == "Yes", ]
vec_species <- subset_df$species 

# Create a list of plots for each species
cols_colours <- c("#00938d","#36528f","#4b0257","#b91a1a")
plots <- lapply(vec_species, function(species) {
  ggplot(df_lm2) +
    aes(x = Risk_group_num, y = df_lm2[[species]]) +
    geom_boxplot(aes(group = Risk_group), outlier.shape = TRUE, fill = cols_colours, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = species, y = "CLR abundance species")  # Optional: Add a title for each plot
})

# Combine the plots using cowplot
combined_plot <- cowplot::plot_grid(plotlist = plots, nrow = 1)
print(combined_plot)


#=====Analyses:pathways=====# ----
#### Import pathway data ----
#Ranko's functions
prepCleanHumann <- function(inPath, dropUnintegrated = T, dropUnmapped = T, dropTaxonSpecific = T,
                            presenceFilter = -1,minRelativeAbundance = -1, rescaleTaxa = T, novogeneIdsClean = T) {
  inDF <- read.table(inPath,sep='\t',header=T,quote ='',comment.char = '')
  #fix pathway ID (these tend to be weird coming out of humann)
  rownames(inDF) <- inDF$X..Pathway
  inDF$X..Pathway <- NULL
  #drop "junk" from humann (unintegrated/unmapped data & taxon-specific pathways)
  if (dropTaxonSpecific) {
    inDF <- inDF[grep('\\|',rownames(inDF),invert = T),]
  }
  if (dropUnintegrated) {
    inDF <- inDF[grep('UNINTEGRATED',rownames(inDF),invert = T),]
  }
  if (dropUnmapped) {
    inDF <- inDF[grep('UNMAPPED',rownames(inDF),invert = T),]
  }
  rownames(inDF)[grep('PWY',rownames(inDF),invert = T)] <- paste0('PWY_',rownames(inDF)[grep('PWY',rownames(inDF),invert = T)])
  inDF <- as.data.frame(t.data.frame(inDF))
  # fix sample IDs, remove duplicates
  inDF$ID <- rownames(inDF)
  if (novogeneIdsClean) {
    print ('NOTE: doing cleaning of Novogene IDs (keeping format <AB>_<CD>_<EFG...>)')
    for (i in c(1:nrow(inDF))) {
      ss <- strsplit(inDF$ID[i],'_')[[1]]
      sss <- paste0(ss[1],'_',ss[2],'_',ss[3])
      inDF$ID[i] <- sss
    }
  }
  if (sum(duplicated(inDF$ID) > 0)) {
    print(paste('WARNING: found ',sum(duplicated(inDF$ID) > 0),'duplicates, dropping them!'))
  }
  inDF <- inDF[!duplicated(inDF$ID),]
  rownames(inDF) <- inDF$ID
  inDF$ID <- NULL
  # make sure columns are actually numbers (otherwise filter dies; NOTE: this should not be necessary, but... )
  for (c in colnames(inDF)) {inDF[[c]] <- as.numeric(inDF[[c]])}
  # clean, rescale and save
  inDFt2 <- filterHumannDF(inDF,presPerc = presenceFilter,minMRelAb = minRelativeAbundance,minMedRelAb = -1,rescale = T,minSum = 1,verbose = T)
  inDFt2 <- inDFt2[,colSums(inDFt2)!=0]
  inDFt2$ID <- row.names(inDFt2)
  print(paste('Done, returning ',nrow(inDFt2),'samples'))
  inDFt2
}
filterHumannDF <- function(inDF,presPerc = 0.05,minMRelAb = 0.001,minMedRelAb=0.0,minSum=90.0, rescale=T,verbose=T,type='MetaCyc') {
  
  nonPWYpwys <- c("ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis",
                  "CHLOROPHYLL-SYN: chlorophyllide a biosynthesis I (aerobic, light-dependent)",
                  "GLYCOLYSIS-E-D: superpathway of glycolysis and Entner-Doudoroff",
                  "GLYCOLYSIS-TCA-GLYOX-BYPASS: superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
                  "GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)",
                  "GLYOXYLATE-BYPASS: glyoxylate cycle",
                  "HEME-BIOSYNTHESIS-II: heme biosynthesis I (aerobic)",
                  "MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS: protein N-glycosylation (eukaryotic, high mannose)",
                  "NAD-BIOSYNTHESIS-II: NAD salvage pathway II",                  
                  "REDCITCYC: TCA cycle VIII (helicobacter)",
                  "TCA-GLYOX-BYPASS: superpathway of glyoxylate bypass and TCA",
                  "TCA: TCA cycle I (prokaryotic)")
  
  colnames(inDF)[colnames(inDF) %in% nonPWYpwys] <- paste0('PWY_',colnames(inDF)[colnames(inDF) %in% nonPWYpwys])
  
  if (type=='MetaCyc') {
    nonPWYdf <- as.data.frame(inDF[,-grep('PWY',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^EC_',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    nonPWYdf <- as.data.frame(inDF[,-grep('RXN',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^PF[01]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^GO',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^K[012]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^K[012]',colnames(inDF))] ])
  }
  colnames(nonPWYdf) <- cnsNonPWYdf
  if (type=='MetaCyc') {
    yesPWYdf <- as.data.frame(inDF[,grep('PWY',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    yesPWYdf <- as.data.frame(inDF[,grep('^EC_',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    yesPWYdf <- as.data.frame(inDF[,grep('RXN',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    yesPWYdf <- as.data.frame(inDF[,grep('^PF[01]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    yesPWYdf <- as.data.frame(inDF[,grep('^GO',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    yesPWYdf <- as.data.frame(inDF[,grep('^K[012]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^K[012]',colnames(inDF))] ])
  }
  
  # replaces NAs with 0s
  for (c in colnames(yesPWYdf)) {
    yesPWYdf[,c][is.na(yesPWYdf[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    nrnZ = as.numeric(sum(yesPWYdf[,c]!=0.0))
    if (nrnZ/as.numeric(nrow(yesPWYdf)) < presPerc) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = mean(yesPWYdf[,c])
    if ( mn < minMRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = median(yesPWYdf[,c])
    if ( mn < minMedRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # do final rescale
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  inDF <- cbind.data.frame(nonPWYdf,yesPWYdf)
  if (verbose) {print ('> DONE')}
  inDF
}

#Clean PWY with functions above
Cleaned_PWY_STEMI <- prepCleanHumann("~/Desktop/results_humann3_pathways_duplicated_removed.tsv",
                                     dropUnintegrated = T,dropUnmapped = T, dropTaxonSpecific = T,
                                     presenceFilter = -1, minRelativeAbundance = -1, rescaleTaxa = T, novogeneIdsClean = T) #keeps 440 pathways

Cleaned_PWY_DAG3 =read.csv("~/Desktop/DAG3_humann3_headersfixed_transposed_cleaned_normalized_v2.csv") #has 606 pathways
row.names(Cleaned_PWY_DAG3) <- Cleaned_PWY_DAG3$ID

#filter participants for analyses
Cleaned_PWY_STEMI %>% filter(ID %in% Meta_merged$sample_ID) -> Cleaned_PWY_STEMI_filtered
Cleaned_PWY_STEMI_filtered$ID <- NULL
Cleaned_PWY_DAG3 %>% filter(ID %in% Meta_merged$sample_ID) -> Cleaned_PWY_DAG3_filtered #two LL particiapants are not in the PWY dataframe ("LL68_E05_6363" "LL79_H09_7454")
Cleaned_PWY_DAG3_filtered$ID <- NULL

#filter out pathways that are zero in filtered participants
Cleaned_PWY_STEMI_filtered <- Cleaned_PWY_STEMI_filtered[,colSums(Cleaned_PWY_STEMI_filtered)>0]
Cleaned_PWY_DAG3_filtered <- Cleaned_PWY_DAG3_filtered[,colSums(Cleaned_PWY_DAG3_filtered)>0]

## Colnames do not overlap, have to make them the same
colnames(Cleaned_PWY_STEMI_filtered) <- gsub("PWY_", "", colnames(Cleaned_PWY_STEMI_filtered))
colnames(Cleaned_PWY_STEMI_filtered) <- str_replace_all(colnames(Cleaned_PWY_STEMI_filtered), "[^[:alnum:]]", ".")
Cleaned_PWY_STEMI_filtered %>% rename("X1CMET2.PWY..N10.formyl.tetrahydrofolate.biosynthesis" = "1CMET2.PWY..N10.formyl.tetrahydrofolate.biosynthesis") -> Cleaned_PWY_STEMI_filtered
Cleaned_PWY_STEMI_filtered %>% rename("X3.HYDROXYPHENYLACETATE.DEGRADATION.PWY..4.hydroxyphenylacetate.degradation" = "3.HYDROXYPHENYLACETATE.DEGRADATION.PWY..4.hydroxyphenylacetate.degradation") -> Cleaned_PWY_STEMI_filtered

#Merge data cohorts
Cleaned_PWY_STEMI_filtered$ID <- rownames(Cleaned_PWY_STEMI_filtered)
Cleaned_PWY_DAG3_filtered$ID <- rownames(Cleaned_PWY_DAG3_filtered)
PWYs_merged <- full_join(Cleaned_PWY_STEMI_filtered, Cleaned_PWY_DAG3_filtered) #keeps 509 pathways in df
rownames(PWYs_merged) = PWYs_merged$ID # making ID the rowname again
PWYs_merged$ID <- NULL

#Removing pathways that are not in both cohorts
PWYs_merged %>% select_if(~ !any(is.na(.))) -> PWYs_merged #keeps 417 pathways

#### CLR transformation (after filtering on prevalence) ----
#First filtering on prevalence (10%)
Prevalence.filter <- function(x, threshold = 10){
  remove_cols <- vector()
  for (i in 1:ncol(x)) {
    cname <- colnames(x)[i]
    if(colSums(x[i] != 0) / nrow(x) *100 < threshold){
      remove_cols <- c(remove_cols, cname) 
    }
  }
  x_new <- x %>% dplyr::select(-remove_cols)
  print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))
  return(x_new)
}
PWYs_filter <- Prevalence.filter(PWYs_merged, threshold = 10)

#Do CLR transformation
my_pseudocount=min(PWYs_filter[PWYs_filter!=0])/2.
PWYs_CLR <- decostand(PWYs_filter, method = "clr", pseudocount = my_pseudocount)

#### Beta diversity calculations and plotting ----
#Calculating beta diversity
vegdist(PWYs_CLR, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- merge(PC, Meta_merged, by = 'row.names')

#Making plots
PCoA_meta$Risk_group <- ordered(PCoA_meta$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI"))

ref <- reformulate("Risk_group","cbind(V1,V2,V3,V4,V5)")
centroids <- aggregate(ref,PCoA_meta,mean) #calculate centroids

my_cols <- c("STEMI"="#b91a1a", "Low"="#00938d", "Intermediate"="#36528f", "High"="#4b0257")

PCoA1_2 <- PCoA_meta %>% ggplot(aes(x=V1, y=V2, color=Risk_group)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_2 <- PCoA1_2 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V2"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V2",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA1_2)

PCoA1_3 <- PCoA_meta %>% ggplot(aes(x=V1, y=V3, color=Risk_group)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_3 <- PCoA1_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V3",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA1_3)

PCoA1_4 <- PCoA_meta %>% ggplot(aes(x=V1, y=V4, color=Risk_group)) + 
  xlab(paste("PCo1 (", var_expl[1],"%)", sep="")) +
  ylab(paste("PCo4 (", var_expl[4],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA1_4 <- PCoA1_4 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V1",y="V4"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V1",y="V4",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA1_4)

PCoA2_3 <- PCoA_meta %>% ggplot(aes(x=V2, y=V3, color=Risk_group)) + 
  xlab(paste("PCo2 (", var_expl[2],"%)", sep="")) +
  ylab(paste("PCo3 (", var_expl[3],"%)", sep="")) +
  geom_point(size=3,alpha=0.4) + theme_bw() + stat_ellipse() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = my_cols)
PCoA2_3 <- PCoA2_3 + geom_point(data=centroids,shape=16,stroke=3,size=4,aes_string(x="V2",y="V3"),alpha=1) +
  geom_point(data=centroids,shape=21,stroke=2,size=4,aes_string(x="V2",y="V3",col="Risk_group"),alpha=0.8, colour="black") 
print(PCoA2_3)

Combined_PCoAs <- ggarrange(PCoA1_2, PCoA1_3, PCoA1_4, PCoA2_3, ncol = 2, nrow = 2)
print(Combined_PCoAs)

#test differences between centroids
wilcoxVarsTouse <- c("V1","V2","V3","V4","V5")
Wilcox_PCoA_Results <- NULL
for (i in wilcoxVarsTouse[1:5]) {
  #create dataframe with only relevant info
  PCoA_meta %>% dplyr::select(c(i, Risk_group)) -> tmDF
  #do wilcoxon
  pairwise.wilcox.test(tmDF[,1], tmDF$Risk_group, p.adjust.method="none") -> wilcox_df
  melt(wilcox_df$p.value) -> wilcoxon_coordinates
  #save info in a dataframe
  Wilcox_all <- data.frame(Coordinate=i,
                           Variable1=wilcoxon_coordinates[,1],
                           Variable2=wilcoxon_coordinates[,2],
                           P_value=wilcoxon_coordinates[,3])
  print(Wilcox_all)
  Wilcox_PCoA_Results <- rbind.data.frame(Wilcox_PCoA_Results,Wilcox_all)
}

#### Adonis (univariate) ----
# load phenotypes
inPhenos <- Meta_merged
# make character variables factors
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)
ordered(inPhenos$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI")) -> inPhenos$Risk_group
str(inPhenos) #check

# remove IDs and other (irrelevant) phenotypes
adonisVarsTouse <- colnames(inPhenos)
toExclude = c("sample_ID", "sex", "age", "sys_bp", "dys_bp", "diabetes", "oral_antidiabetic", "insulin", "hdl", "cholesterol_total", "ldl",
              "smoking_status", "at_antagonist", "beta_block", "ca_antagonist", "ace_inhibitor","diuretics_low_ceiling")
adonisVarsTouse <- adonisVarsTouse[!adonisVarsTouse %in% toExclude]

adonisResults <- NULL
permNR <- 1000
inMB <- PWYs_CLR
inMB <- tibble::rownames_to_column(inMB, "sample_ID")

for (i in adonisVarsTouse[1:18]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"sample_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="sample_ID",y=inMB,by.y="sample_ID")
  rownames(allDF) <- allDF$sample_ID
  allDF$sample_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis(inBC ~ av,permutations=permNR)
    aov_table <- ad$aov.tab
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov_table[1,1],
                         SumsOfSqs=aov_table[1,2],
                         MeanSqs=aov_table[1,3],
                         FModel=aov_table[1,4],
                         R2=aov_table[1,5],
                         pval=aov_table[1,6],
                         FDR.BH=NA,
                         Significant=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    write.table(adonisResults,paste0("~/Documents/MDPhD/Results_STEMI_DAG3/adonis_results_table.csv"),sep=",",row.names=F)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$Var <- NULL
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant="No"
adonisResults$Significant[adonisResults$FDR.BH<0.05]="Yes"

#plot for PWYs
adonis_PYWs <- ggplot(adonisResults, aes(reorder(row.names(adonisResults), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Pathways") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
plot(adonis_PYWs)

#### Adonis (multivariate) ----
All_df <- merge(Meta_merged, PWYs_CLR, by = 'row.names')
ad_meta <- All_df[c(1:37)]
ad_pwy <- All_df[c(38:408)]
ad_pwy_dm <- vegdist(ad_pwy,method = "euclidean")
ad2 <- adonis2(formula = ad_pwy_dm ~ proton_pump_inhib + Risk_group + reads, data = ad_meta, permutations = 1000, by = "margin")

model.frame(ad2) -> adonis_mv
adonis_mv$Significant="No"
adonis_mv$Significant[adonis_mv$`Pr(>F)`<0.05]="Yes"

adonis_PWYs_multi <- ggplot(adonis_mv, aes(reorder(row.names(adonis_mv), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariate: Pathways") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
plot(adonis_PWYs_multi)

#### Emmeans per contrast FDR ----
df_lm2 <- Meta_merged[c("sample_ID", "proton_pump_inhib", "reads", "Risk_group")]
df_lm2[sapply(df_lm2, is.character)] <- lapply(df_lm2[sapply(df_lm2, is.character)], as.factor)
df_lm2$sample_ID <- as.character(df_lm2$sample_ID)
df_lm2$Risk_group <- factor(df_lm2$Risk_group, ordered = FALSE)
df_lm2$Risk_group <- relevel(df_lm2$Risk_group, ref = "Low") 
levels(df_lm2$Risk_group)
PWYs_CLR$sample_ID <- rownames(PWYs_CLR)

df_lm2 <- left_join(df_lm2, PWYs_CLR, by="sample_ID")

myContr1 <-list(Low_Intermediate=c(-1,0,1,0),
                Low_High=c(-1,1,0,0),
                Low_STEMI=c(-1,0,0,1),
                Intermediate_High=c(0,1,-1,0),
                Intermediate_STEMI=c(0,0,-1,1),
                High_STEMI=c(0,-1,0,1),
                LIH_STEMI=c(-(1/3),-(1/3),-(1/3),1))

lm_contrasts <- vector("list", ncol(PWYs_CLR))
names(lm_contrasts) <- colnames(PWYs_CLR)
lm_models <- vector("list", ncol(PWYs_CLR))
names(lm_models) <- colnames(PWYs_CLR)

for(i in colnames(df_lm2)[5:375]) {
  #tryCatch({
  print(i)
  m_pathways <- lm(df_lm2[,i] ~ Risk_group + proton_pump_inhib + reads, 
                  data=df_lm2[,c(i,"Risk_group", "proton_pump_inhib", "reads")])
  lm_models[[i]] <- m_pathways
  mm_pathways <- emmeans::emmeans(m_pathways,"Risk_group")
  lm_contrasts[[i]] <- as.data.frame(emmeans::contrast(mm_pathways, myContr1)) } #transform the emmGrid object into a data frame, which will be recognized as a vector in tibble

#Get all results in one dataframe
lm_raw_results <- lm_contrasts %>% bind_rows(.id="pathway_id") 

#Calculate FDR per contrast
lm_contrasts_sigFDR <- vector("list",7)
names(lm_contrasts_sigFDR) <- c("Low_Intermediate", "Low_High","Low_STEMI","Intermediate_High","Intermediate_STEMI","High_STEMI","LIH_STEMI")
lm_contrasts_allFDR <- vector("list",7)
names(lm_contrasts_allFDR) <- c("Low_Intermediate", "Low_High","Low_STEMI","Intermediate_High","Intermediate_STEMI","High_STEMI","LIH_STEMI")

for(c in c("Low_Intermediate", "Low_High","Low_STEMI","Intermediate_High","Intermediate_STEMI","High_STEMI","LIH_STEMI")) {
  df <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    filter(p_adj < 0.05) #%>%
    #mutate(pathway_id=paste0(str_split_fixed(pathway_id,"\\.",7)[,7]),
           # sign=ifelse(estimate>0, "positive", "negative"),
          # contrast=factor(contrast))
  lm_contrasts_sigFDR[[c]] <- df
  
  df1 <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH"))# %>%
   # mutate(pathway_id=paste0(str_split_fixed(pathway_id,"\\.",7)[,7]),
           # sign=ifelse(estimate>0, "positive", "negative"),
          # contrast=factor(contrast))
  lm_contrasts_allFDR[[c]] <- df1
}

lm_contrasts_sigFDR_PWY <- lm_contrasts_sigFDR %>% bind_rows() #significant dataframe
lm_contrasts_allFDR_PWY <- lm_contrasts_allFDR %>% bind_rows() #all results

write.xlsx(lm_contrasts_allFDR_PWY,"~/Documents/MDPhD/Results_STEMI_DAG3/DDA_pathways.xlsx")


#Looking deeper into the pwys
unique(lm_contrasts_sigFDR_PWY$pathway_id)

library(readxl)
Metacyc_pathways <- read_excel("Documents/MDPhD/Metacyc pathways.xlsx")
Metacyc_pathways %>% filter(`MetaCyc ID` %in% lm_contrasts_sigFDR_PWY$pathway_id) -> PWY_info
PWY_info %>% dplyr::select(c(1,8:18)) -> PWY_info
rename(PWY_info, pathway_id = `MetaCyc ID`) -> PWY_info
full_join(lm_contrasts_sigFDR_PWY, PWY_info) -> PWY_sig_info

unique(PWY_sig_info$pathway_id)

ggplot(df_lm2, aes(x = Risk_group, y = PWY.622..starch.biosynthesis)) +
  geom_boxplot()


#### Heatmap plotting (pathway) ----
complete(lm_contrasts_sigFDR_PWY, pathway_id, contrast) -> df_heatmap_pathways
rng_estimates <- range(df_heatmap_pathways[!is.na(df_heatmap_pathways$estimate),]$estimate)

full_join(df_heatmap_pathways, PWY_info, by="pathway_id") -> df_heatmap_info

ggplot(df_boxplot, aes(x=Risk_group, y=df_boxplot[,i])) +
  geom_boxplot(aes(fill=Risk_group, alpha=0.6))

heatmap_pathways <- df_heatmap_info %>% 
  ggplot(aes(y=pathway_id,x=contrast,fill=estimate)) +
  coord_fixed()+ # make the cell contain square
  geom_tile(colour="black") +
  scale_fill_gradient2(name="Effect size", low = "#1111AA", mid = "white", high = "#AA1111",midpoint = 0, na.value = "grey90",limits=c(floor(rng_estimates[1]), ceiling(rng_estimates[2]))) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.text.y = element_text(size=8, face="italic", colour = "black"),
        axis.text.x = element_text(size=6, face="bold", colour = "black", hjust= 1, angle=45),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.position = "top",
        legend.text = element_text(size=6, face="bold", colour = "black"),
        legend.title =element_text(size=6, face="bold", colour = "black"),
        legend.margin=margin(grid::unit(0, "cm")),
        legend.key.height=grid::unit(0.2, "cm"),
        legend.key.width=grid::unit(0.8, "cm")) +
  labs(x=NULL, y=NULL)
print(heatmap_pathways)


#play around
# Define the data frame with the relevant columns
df <- df_heatmap_info[, c("pathway_id", "contrast", "estimate", "Category")]
# Reshape the data frame to create a matrix
mat <- dcast(df, pathway_id ~ contrast, value.var = "estimate")
rownames(mat) <- mat$pathway_id
mat$pathway_id <- NULL
mat[is.na(mat)] <- 0

ComplexHeatmap::Heatmap(mat)
str(mat)



# Define the row annotations
row_anno <- data.frame(Category = df$Category[match(row.names(mat), df$pathway_id)], row.names = row.names(mat))
row_anno <- data.frame(Category = ifelse(is.na(row_anno$Category), "no data", row_anno$Category))
row_anno$Category <- as.factor(row_anno$Category)

# Create the heatmap with row annotations
pheatmap(mat, annotation_row = row_anno, cluster_rows = FALSE)

library(ComplexHeatmap)

#====Other====#----
#### Table 1 ----
Meta_merged$Risk_group <- ordered(Meta_merged$Risk_group, levels = c("Low", "Intermediate", "High", "STEMI"))

table1(~ age + sex + bmi + dys_bp + sys_bp + heart_rate + 
         diabetes + hypercholesterolemia + hypertension + smoking_status + arterydisease + CVA + 
         leukocytes + creatinine + triglyceride + hdl + ldl + cholesterol_total +
         statin + proton_pump_inhib + beta_block + at_antagonist + ca_antagonist + oral_antidiabetic +
         diuretics_low_ceiling + ace_inhibitor + insulin + antiarrythmica + inhalation + lipid_lowering +
         P2Y12_inhibitors + anti_trombotics + reads + Framingham| Risk_group, data=Meta_merged)

#Test to see what variables differ significantly between groups
Meta_merged %>% dplyr::select(Risk_group, everything()) -> Meta_merged #make risk groups the first column
Meta_merged %>% dplyr::select(c(1:34, 36)) -> Metadata_kruskal

Kruskal_df <- NULL
for(i in colnames(Metadata_kruskal[3:35])) {
  print(i)
  mF <- formula(paste(i, " ~ Risk_group"))
  my_test <- kruskal.test(mF, data = Metadata_kruskal)
  # accumulate results
  oneRow <- data.frame(Variable=i,
                       Statistic=my_test$statistic,
                       Parameter=my_test$parameter,
                       Pval=my_test$p.value,
                       FDR.BH=NA,
                       Significant=NA)
  print(oneRow)
  Kruskal_df <- rbind.data.frame(Kruskal_df,oneRow) 
}

Kruskal_df$FDR.BH=p.adjust(Kruskal_df$Pval, method = "BH")
Kruskal_df$Significant="No"
Kruskal_df$Significant[Kruskal_df$FDR.BH<0.05]="Yes"

write.xlsx(Kruskal_df,"~/Documents/MDPhD/Results_STEMI_DAG3/Kruskal_Wallis.xlsx")

