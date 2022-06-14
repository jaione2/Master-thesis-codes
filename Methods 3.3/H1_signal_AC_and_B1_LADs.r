#--------------------------------------------------------------------------#
# Plot the abundance of CHiPSeq H1 variants in B1 LADs identified by SICER #
#--------------------------------------------------------------------------#
library (ggplot2)
library (reshape2)

#Set working directory 
#Upload data:
data <- as.data.frame(read.table("abundances_H1_B1LADs_tab.bed", header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#Specify the LADs type in which signal has been mapped 
lad<-"B1"

#Add names to the columns

colnames(data) <- c("Chr", "Start", "End", "H1.0Dox", "H1.0NoDox","H1.2Dox", "H1.2NoDox","H1.4Dox", "H1.4NoDox", 
                    "H1.5Dox", "H1.5NoDox","H1.XDox", "H1.XNoDox")

##################
###PLOT THE H1s###
##################
#Reconstruct the data frame: For this, we will take separately the data regardign to each histone
#and add a column with the name of the variant and melt the abudnance data that will be dependent of the "variable" column
#H1.0
dataH1.0<-data[,c(4:5)]
dataH1.0["Histone_variant"]="H1.0"
data.m <- melt(dataH1.0,id.vars="Histone_variant", measure.vars=c('H1.0Dox','H1.0NoDox'))

#H1.2
dataH1.2<-data[,c(6:7)]
dataH1.2["Histone_variant"]="H1.2"
data.m2 <- melt(dataH1.2,id.vars="Histone_variant", measure.vars=c('H1.2Dox','H1.2NoDox'))

#h1.4
dataH1.4<-data[,c(8:9)]
dataH1.4["Histone_variant"]="H1.4"
data.m3 <- melt(dataH1.4,id.vars="Histone_variant", measure.vars=c('H1.4Dox','H1.4NoDox'))
#H1.5
dataH1.5<-data[,c(10:11)]
dataH1.5["Histone_variant"]="H1.5"
data.m4 <- melt(dataH1.5,id.vars="Histone_variant", measure.vars=c('H1.5Dox','H1.5NoDox'))
#H1X
dataH1X<-data[,c(12:13)]
dataH1X["Histone_variant"]="H1X"
data.m5 <- melt(dataH1X,id.vars="Histone_variant", measure.vars=c('H1.XDox','H1.XNoDox'))

#Bind all the information
histone_data<-rbind(data.m,data.m2,data.m3,data.m4,data.m5)

#Translate the Dox/NoDox information from the "variable" column into a new binomial column called "treatment"
histone_data$treatment<-0
histone_data$treatment[histone_data$variable=="H1.0Dox"]<-"2Dox"
histone_data$treatment[histone_data$variable=="H1.2Dox"]<-"2Dox"
histone_data$treatment[histone_data$variable=="H1.4Dox"]<-"2Dox"
histone_data$treatment[histone_data$variable=="H1.5Dox"]<-"2Dox"
histone_data$treatment[histone_data$variable=="H1.XDox"]<-"2Dox"

histone_data$treatment[histone_data$variable=="H1.0NoDox"]<-"1NoDox"
histone_data$treatment[histone_data$variable=="H1.2NoDox"]<-"1NoDox"
histone_data$treatment[histone_data$variable=="H1.4NoDox"]<-"1NoDox"
histone_data$treatment[histone_data$variable=="H1.5NoDox"]<-"1NoDox"
histone_data$treatment[histone_data$variable=="H1.XNoDox"]<-"1NoDox"
table(histone_data$treatment)

#Plot the H1s signal abundance comparing NoDox/Dox treatment
lower<- -0.01  
upper<- 0.01

title<-paste("ChIPsSeq signal of H1 variants in ", lad, " LADs found with EDD", sep="")

ggplot(histone_data)+
  geom_boxplot(mapping=aes(x=Histone_variant,y=value, fill=treatment),outlier.shape = NA)+
  scale_y_continuous(limits = c(lower, upper))+
  theme_bw()+
  labs(title=title,x="",y="ChIPSeq signal")+
  theme(plot.title = element_text(hjust=0.5, size=14),axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 12,color="black"),legend.text = element_text(size = 12))+
  scale_x_discrete(labels=c("H1.0","H1.2","H1.4","H1.5","H1X"))+
  scale_fill_manual(values = c("#E69F00","#56B4E9"),label=c("No Dox","Dox"),name="Treatment")+
  theme(legend.position = "right",legend.background = element_rect(fill = "white", color = "black"))
