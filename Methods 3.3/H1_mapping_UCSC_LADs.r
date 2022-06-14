#Plot the abundance of CHiPSeq reads mapping in the G bands 
library (ggplot2)
library(reshape2)

#Set working directory 
#We must upload 5 files 

data_ilad<-as.data.frame(read.table("signals_H1_ILAD.bed", header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
data_ucsc<-as.data.frame(read.table("signals_H1_ucsc.bed", header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
data_random<-as.data.frame(read.table("signals_H1_random.bed", header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#As it has been random, there are artificial chr in the analyisis and they will 
#be removed in this step 
data_random <- data_random[data_random$V4 != ".",]
nrow(data_random)
library(qpcR)
#Prepare the separate dataframe for each H1 

#Add colnames to the dataframes
colnames(data_ilad)<-c("Chr", "Start", "End","H1.0","H1.2","H1.4","H1.5","H1X")
colnames(data_ucsc)<-c("Chr", "Start", "End","H1.0","H1.2","H1.4","H1.5","H1X")
colnames(data_random)<-c("Chr", "Start", "End","H1.0","H1.2","H1.4","H1.5","H1X")

#H1.0
data.0<-qpcR:::cbind.na(data_ilad[,4],data_ucsc[,4],data_random[,4])
data.0<-as.data.frame(data.0)
colnames(data.0)<-c("ILAD","UCSC","RANDOM")
data.0["Histone"]="H1.0"
class(data.0)
data_n0<- melt(data.0,id.vars="Histone", measure.vars=c("ILAD","UCSC","RANDOM"))

#H1.2
data.1<-qpcR:::cbind.na(data_ilad[,5],data_ucsc[,5],data_random[,5])
data.1<-as.data.frame(data.1)
colnames(data.1)<-c("ILAD","UCSC","RANDOM")
data.1["Histone"]="H1.2"
class(data.1)
data_n2<- melt(data.1,id.vars="Histone", measure.vars=c("ILAD","UCSC","RANDOM"))

#H1.4
data.2<-qpcR:::cbind.na(data_ilad[,6],data_ucsc[,6],data_random[,6])
data.2<-as.data.frame(data.2)
colnames(data.2)<-c("ILAD","UCSC","RANDOM")
data.2["Histone"]="H1.4"
class(data.2)
data_n4<- melt(data.2,id.vars="Histone", measure.vars=c("ILAD","UCSC","RANDOM"))

#H1.5
data.3<-qpcR:::cbind.na(data_ilad[,7],data_ucsc[,7],data_random[,7])
data.3<-as.data.frame(data.3)
colnames(data.3)<-c("ILAD","UCSC","RANDOM")
data.3["Histone"]="H1.5"
class(data.3)
data_n5<- melt(data.3,id.vars="Histone", measure.vars=c("ILAD","UCSC","RANDOM"))

#H1X 
data.4<-qpcR:::cbind.na(data_ilad[,8],data_ucsc[,8],data_random[,8])
data.4<-as.data.frame(data.4)
colnames(data.4)<-c("ILAD","UCSC","RANDOM")
data.4["Histone"]="H1X"
class(data.4)
data_nx<- melt(data.4,id.vars="Histone", measure.vars=c("ILAD","UCSC","RANDOM"))

histone_data<-rbind(data_n0,data_n2,data_n4,data_n5,data_nx)

histone_data$treatment<-0
histone_data$treatment[histone_data$variable=="ILAD"]<-"ILAD"
histone_data$treatment[histone_data$variable=="UCSC"]<-"LAD"
histone_data$treatment[histone_data$variable=="RANDOM"]<-"random"

nrow(histone_data)
nrow(na.omit(histone_data))

histone_data<-na.omit(histone_data)
nrow(histone_data)
class(histone_data$value)
histone_data$value<-as.numeric(histone_data$value)

title<-paste("ChIPSeq signal of T47D NoDox H1s in UCSC LADs",sep="")

ggplot(histone_data)+
  geom_boxplot(mapping=aes(x=Histone,y=value,fill=treatment),outlier.shape=NA)+
  labs(title=title,x="",y="CHiPSeq signal")+
  ylim(-0.02,0.03)+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5, size=14),axis.text.y=element_text(margin=margin(t=0, r=15, b=0, l=0), hjust=0.5, size=14),axis.text.x=element_text(size=14), strip.text = element_text(size=13),legend.title=element_text(size=10), legend.text=element_text(size=10))+
  scale_x_discrete(labels=c("H1.0","H1.2","H1.4","H1.5","H1X"))+
  scale_fill_discrete(name="",label=c("Inter LADs","LADs","Random LADs"))+
  theme(legend.position="right")