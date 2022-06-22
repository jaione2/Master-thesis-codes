##H1 variants abundance in SKNSH cells  at chromatin colors classified by CHROMHMM
#ChIP-seq data: H1 variants in SKNSH cells from ChIPSeq 9b and BGI
#Chromatin segments  download from UCSC genome browser --> ChromHMM HeLaS3

#directory 

setwd("/home/ibmb-ajv-03/Bioinformatics_Jaione/Chromatin states heatmap/ChIPSeq_T47DNoDox_ChromColors/")

#Load Files data. ChIPseq5b (CRG) data mapped (Bedtools map) at ChromHMM states
#dim 1421299x5
H1.0<-read.table("H1.0NoDox_T47D6c2_mapped_at_ChromHMMHeLaS3.bed",stringsAsFactors = FALSE)
H1.2 <- read.table("H1.2NoDox_Oct_T47D6c2_mapped_at_ChromHMMHeLaS3.bed", stringsAsFactors = FALSE)
H1.4 <- read.table("H1.4NoDox_L1L2Merged_T47D6c_mapped_at_ChromHMMHeLaS3.bed", stringsAsFactors = FALSE)
H1.5 <- read.table("H1.5NoDox_T47D6c2_mapped_at_ChromHMMHeLaS3.bed", stringsAsFactors = FALSE)
H1X <- read.table("H1XNoDox_Oct_T47D6c2_mapped_at_ChromHMMHeLaS3.bed", stringsAsFactors = FALSE)


#####Define histones to plot
histones_vector <- c("H1.0","H1.2", "H1.4", "H1.5", "H1X")
histones_list <- list(H1.0, H1.2, H1.4, H1.5, H1X)
names(histones_list) <- c(histones_vector)
#########


#Merge data to create a single data frame
H1variants_ChromHMM <- histones_list[[1]][,1:4]

for (histone in histones_vector){
  H1variants_ChromHMM <- cbind(H1variants_ChromHMM, histones_list[[histone]][,5])
}

colnames(H1variants_ChromHMM) <- c("chr", "start","end","ChromState",histones_vector)


#Agruppate ChromStates to 10 groups (as defined per ENCODE) iwth complete name

H1variants_ChromHMM$ChromGroup <- ifelse(H1variants_ChromHMM$ChromState=="Tss"|H1variants_ChromHMM$ChromState=="TssF","Active_Promoter",
                                         ifelse(H1variants_ChromHMM$ChromState=="PromF","Promoter_flanking",
                                                ifelse(H1variants_ChromHMM$ChromState=="PromP","Inactive_Promoter",
                                                       ifelse(H1variants_ChromHMM$ChromState=="Enh"|H1variants_ChromHMM$ChromState=="EnhF","Candidate_strong_enhancer",
                                                              ifelse(H1variants_ChromHMM$ChromState=="EnhWF"|H1variants_ChromHMM$ChromState=="EnhW"|H1variants_ChromHMM$ChromState=="DnaseU"|H1variants_ChromHMM$ChromState=="DnaseD"|H1variants_ChromHMM$ChromState=="FaireW","Candidate_weak_enhancer",
                                                                     ifelse(H1variants_ChromHMM$ChromState=="CtcfO"|H1variants_ChromHMM$ChromState=="Ctcf","Candidate_Insulator",
                                                                            ifelse(H1variants_ChromHMM$ChromState=="Gen5'"|H1variants_ChromHMM$ChromState=="Elon"|H1variants_ChromHMM$ChromState=="ElonW"|H1variants_ChromHMM$ChromState=="Gen3'"|H1variants_ChromHMM$ChromState=="Pol2"|H1variants_ChromHMM$ChromState=="H4K20","Transcription-associated",
                                                                                   ifelse(H1variants_ChromHMM$ChromState=="Low","Low",
                                                                                          ifelse(H1variants_ChromHMM$ChromState=="ReprD"|H1variants_ChromHMM$ChromState=="Repr"|H1variants_ChromHMM$ChromState=="ReprW","Polycomb-repressed",ifelse(H1variants_ChromHMM$ChromState=="Quies"|H1variants_ChromHMM$ChromState=="Art","Heterochr","None")
                                                                                          )))))))))     


#Save RFile
save(H1variants_ChromHMM,file="H1abundance_t47d_nodox_ChromStatesChromHMM_HeLaS3.RData")


###################################################################################
###################################################################################

#10 chromatin groups as factors
groups <- factor(H1variants_ChromHMM$ChromGroup,
                 levels=c("Active_Promoter", "Promoter_flanking", "Inactive_Promoter", "Candidate_strong_enhancer","Candidate_weak_enhancer",
                          "Candidate_Insulator", "Transcription-associated", "Low", "Polycomb-repressed", "Heterochr"))

colors_Boxplot <- c("Active_Promoter"="red1","Promoter_flanking"="tomato","Inactive_Promoter"="violetred3","Candidate_strong_enhancer"="darkorange",
                    "Candidate_weak_enhancer"="yellow1","Candidate_Insulator"="dodgerblue","Transcription-associated"="seagreen4",
                    "Low"="seagreen1", "Polycomb-repressed"="gray20", "Heterochr"="gray")

###################################################################################
###################################################################################


##################################################
######Delete outliers 


summary(H1variants_ChromHMM)  #me hago idea de los Max y Min para compararlos una vez quitado los extremos

#defino los limites: 1st and 3rd Qu.
UpLim <- c(); DownLim <- c()
for (histone in histones_vector){
  DownLim <- c(DownLim, summary(H1variants_ChromHMM[,histone])[[2]]-0.04)
  UpLim <- c(UpLim, summary(H1variants_ChromHMM[,histone])[[5]]+ 0.04)
}
UpLim <- mean(UpLim); DownLim <- mean(DownLim)

#Me quedo con todo lo que esta por debajo del UpLim y por encima del DownLim: POR FILAS COMPLETAS
dim(H1variants_ChromHMM) 
for (histone in histones_vector){
  H1variants_ChromHMM <- H1variants_ChromHMM[(H1variants_ChromHMM[,histone]<UpLim),]
  H1variants_ChromHMM <- H1variants_ChromHMM[(H1variants_ChromHMM[,histone]>DownLim),]
}
# Checking how many rows have been removed
dim(H1variants_ChromHMM) #1245314       

#How many rows have been removed (in %)??
100-((nrow(H1variants_ChromHMM)*100)/1421299)  #12.38%


######################################################
#Data frame with mean values (ChIP-seq abundance and length) per Chrom group (10 groups)
library(dplyr)
Medians_abundance <- H1variants_ChromHMM %>%
  group_by(ChromGroup) %>%
  summarise_at(vars(-chr,-start,-end,-ChromState), list(median=median))

Medians_abundance$ChromGroup_Number <- ifelse(Medians_abundance$ChromGroup=="Active_Promoter", 1,
                                              ifelse(Medians_abundance$ChromGroup=="Promoter_flanking", 2,
                                                     ifelse(Medians_abundance$ChromGroup=="Inactive_Promoter", 3,
                                                            ifelse(Medians_abundance$ChromGroup=="Candidate_strong_enhancer", 4,
                                                                   ifelse(Medians_abundance$ChromGroup=="Candidate_weak_enhancer", 5,
                                                                          ifelse(Medians_abundance$ChromGroup=="Candidate_Insulator", 6,
                                                                                 ifelse(Medians_abundance$ChromGroup=="Transcription-associated", 7,
                                                                                        ifelse(Medians_abundance$ChromGroup=="Low", 8,
                                                                                               ifelse(Medians_abundance$ChromGroup=="Polycomb-repressed", 9,
                                                                                                      ifelse(Medians_abundance$ChromGroup=="Heterochr", 10,"None"))))))))))

Medians_abundance$ChromGroup_Number <- as.numeric(Medians_abundance$ChromGroup_Number)
Medians_abundance <- Medians_abundance[order(Medians_abundance$ChromGroup_Number),]
colnames(Medians_abundance) <- c("ChromGroup", histones_vector, "Length_median","ChromGroup_Number")



######################################################

##Heatmap

# First, pheatmap only takes the numeric matrix object as input.
# So, we need to transfer the numeric part of the data frame to a matrix

numericaldata <- as.matrix(Medians_abundance[,c(2:(ncol(Medians_abundance)-1))])


# An independent data frame with annotations to the rows or columns of the heatmap matrix
# the row names of the annotation data frame have to match the row names or column names 
# of the heatmap matrix depending on your annotation target.

ChromGroup = data.frame("ChromGroup" = Medians_abundance$ChromGroup)

rownames(ChromGroup) <- Medians_abundance$ChromGroup 

rownames(numericaldata) = rownames(ChromGroup)


# List with colors for each annotation.
colors=list(ChromGroup= c("Active_Promoter"="red1","Promoter_flanking"="tomato","Inactive_Promoter"="violetred3","Candidate_strong_enhancer"="darkorange",
                          "Candidate_weak_enhancer"="yellow1","Candidate_Insulator"="dodgerblue","Transcription-associated"="seagreen4",
                          "Low"="seagreen1", "Polycomb-repressed"="gray20", "Heterochr"="gray"))

#Heatmap
library(pheatmap)
library(seriation)
library(dendextend)

#Not scaled
phtmap <- pheatmap(numericaldata, annotation_row = ChromGroup, annotation_colors = colors[1], 
                   cutree_cols = 2,
                   show_rownames = F)
#Scaled
phtmap <- pheatmap(numericaldata, scale="column",annotation_row = ChromGroup, annotation_colors = colors[1], 
                   cutree_cols = 2,
                   show_rownames = F)
col_dend <- phtmap[[2]]
col_dend<-rotate(col_dend, order =c("H1.2", "H1.0", "H1.5", "H1.4", "H1X"))
phtmap
#################
pdf("HeatmapNOScaled_t47d_nodox_median_at_ChromStates.pdf")
pheatmap(numericaldata, cluster_cols=as.hclust(col_dend),  annotation_row = ChromGroup, annotation_colors = colors[1], 
         cutree_cols = 2,cutree_rows = 1,
         show_rownames = F)
dev.off()         
#################


####################################################################################################
########################Select 1000 random fragments of each chromatin state########################
####################################################################################################

#function get 1000 random selection (or whatever number of random selection you indicate). 
#Max around 5000 because there are aprox5000 Inactive promoters (after deleting outliers)
get1000Random <- function(c, h = H1variants_ChromHMM) {
  set.seed(123)
  x <- filter(h, ChromGroup==c)
  x_1000 <- sample(1:nrow(x), 1000)
  x <- x[x_1000,]
}


chrom_states <- c("Active_Promoter", "Promoter_flanking", "Inactive_Promoter", "Candidate_strong_enhancer",
                  "Candidate_weak_enhancer","Candidate_Insulator","Transcription-associated",
                  "Low", "Polycomb-repressed", "Heterochr")


#Create a list in which element is a dataframe containing 1000 random segments of a given chromatin state
random1000List <- list()
for (state in chrom_states) {
  random1000List[[state]] <- get1000Random(state)
}


#Hacer 50 grupos de 20 cada uno (50 si has cogido 1000 random; 250 si has cogido 5000).
#Result: 50random segments de cada chromatin state, cada uno viene de hacer la mediana de 20 random segments.
#Total: 500 valores (50x10chrom states)
juntarMedianasGrupo <- function(grupo, nombre){
  set.seed(125)
  grupo_20_50 <- replicate(50,sample(1:nrow(grupo), 20))
  
  medianasGrupos <- function(d, columnas = 5:(ncol(d)-1)) {
    apply(d[columnas], 2, median)
  }
  
  myList <- list()
  
  for (num in 1:dim(grupo_20_50)[2]) {
    myList[[num]] <- medianasGrupos(grupo[grupo_20_50[,num],])
  }
  
  resp <- as.data.frame(do.call("rbind", myList))
  resp$ChromGroup <- nombre
  return(resp)
}

myCollection <- list()

for (state in random1000List) {
  nombre <- state$ChromGroup[1]
  myCollection[[nombre]] <- juntarMedianasGrupo(state, nombre)
}

result <- do.call("rbind", myCollection)


####################Heatmap 1000 random segments (50x20)
numericaldata_random <- as.matrix(result[,c(1:(ncol(result)-1))])

ChromGroup = data.frame("ChromGroup" = result$ChromGroup)
rownames(ChromGroup) <- rownames(result) 
rownames(numericaldata_random) = rownames(ChromGroup)

colors=list(ChromGroup= c("Active_Promoter"="red1","Promoter_flanking"="tomato","Inactive_Promoter"="violetred3","Candidate_strong_enhancer"="darkorange",
                          "Candidate_weak_enhancer"="yellow1","Candidate_Insulator"="dodgerblue","Transcription-associated"="seagreen4",
                          "Low"="seagreen1", "Polycomb-repressed"="gray20", "Heterochr"="gray"))

#Not scaled
phtmap <- pheatmap(numericaldata_random, annotation_row = ChromGroup, annotation_colors = colors[1],
                   cutree_cols = 2,cutree_rows = 1,
                   show_rownames = F)
#Scaled
phtmap <- pheatmap(numericaldata_random, scale="column",annotation_row = ChromGroup, annotation_colors = colors[1],
                   cutree_cols = 2,cutree_rows = 1,
                   show_rownames = F)
col_dend <- phtmap[[2]]
col_dend <- rotate(col_dend, order =c("H1.5", "H1.5_BGI", "H1.0_BGI","H1.2", "H1.4", "H1.2_BGI", "H1.4_BGI", "H1X_BGI", "H1X"))

################################################################################

###Save Heatmap plot
pdf("Heatmap_NotSCALED_SKNSH_1000randomsegments_median50x20groups_cutcols2_cutrows1.pdf")
pheatmap(numericaldata_random, cluster_cols=as.hclust(col_dend),  annotation_row = ChromGroup, annotation_colors = colors[1],
         cutree_cols = 2,cutree_rows = 1,
         show_rownames = F)
dev.off() 

################################################################################
#####   THE END    #####
################################################################################