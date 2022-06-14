#----------------------------------------------------------#
# Correlation visualization between H1 variants and lamins #
#----------------------------------------------------------#

# In this analysis the aim of is to see the correlation of between the abundance of the H1 variants signals vs. the signal of 
# lamins. The input file is coming from the ChIPSeq_signal_100kb_cor.test.sh. In this case as we had 10 files for the 5
# histones (H1.0,H1.2, H1.4, H1.5 and H1X) in both NoDox and Dox conditions + 4 files for lamins (AC, B1, AC_merged,
# B1_merged). So, the output files have 3 common columns for the coordinates of the 100kb windows (chr, Start, End) + 14
# columns (one for each signal file). 

#------------------------------------------------------------------------------------------------------------------------------#

#Set working directory 
setwd("/home/ibmb-ajv-03/Bioinformatics_Jaione/T47D with ssh (NoDox_Dox)/Plot_H1s_Lamins/")
#Upload the data:
data<-as.data.frame(read.table("abundances_H1_Lamins_tab.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#Change the name of the columns 
colnames(data)
colnames(data) <- c("Chr", "Start", "End", "H1.0Dox", "H1.0NoDox","H1.2Dox", "H1.2NoDox","H1.4Dox", "H1.4NoDox", 
                    "H1.5Dox", "H1.5NoDox","H1.XDox", "H1.XNoDox","LaminB1_merged","Lamin_B1", "LaminAC_merged", "Lamin_AC")
colnames(data)
summary(data)

#SPEARMAN CORRELATION 

# Generate 4 empty vectors for storing data and open 2 loops. The first loop iterates over the H1s and the second over the
# lamins. Store the value obtained from the test in the vectors. 

pvalues<-c()
estimates<-c()
names1<-c()
names2<-c()
for (j in 4:13) { #Iterates over the H1s
  for (k in 14:17){ #Iterate over the lamins
    cor.results<-cor.test(data[,j],data[,k],method="spearman")
    pvalue <- cor.results$p.value
    estimate <- cor.results$estimate
    pvalues<-c(pvalues, pvalue)
    estimates<-c(estimates, estimate)
    names1<-c(names1, colnames(data[j]))
    names2<-c(names2, colnames(data[k]))
  }
}

# Generate a data frame with all the information from the pearson correlation.As the data frame does not allow to repeat the
# same name in rownames, we have to transpose the df.
correlations<-data.frame(names2,estimates,pvalues)
correlations<-as.data.frame(t(correlations))
rownames(correlations)<-c("Lamin type","Estimate","P value")
colnames(correlations)<-names1

# SCATTERPLOT
# Generate 40 plots, 1 per each histone sample vs each lamin sample. The images are downloaded into the file in PNG format. 
library(ggplot2)

for (z in 1:length(names1)){
  histone<-names1[z]
  lamin<-names2[z]
  p.value<-pvalues[z]
  estim<-estimates[z]
  image_name<-paste(histone,"_VS_", lamin,sep="")
  png(image_name, width = 256, height = 192)
  title1<-paste(histone," vs. ",lamin, "\n P value= ",p.value, "\n corr. coef.= ", estim)
  Graph=ggplot(data,aes_string(x=histone,y=lamin)) ### !!! aes_string, si no no pilla
  Graph + geom_point(size=2,shape=21,color="#002344",fill="#FECB00")+
    geom_smooth(method="lm",se=0,size=0.5, color="red")+
    labs(title=title1,x=histone, y=lamin)+
    theme(plot.title = element_text(hjust=0.5, size=12),axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16), axis.text=element_text(size=12))
  dev.off()
}



