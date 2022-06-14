work_dir<-"/home/ibmb-ajv-03/Bioinformatics_Jaione/T47D Dox/UCSC_intersection/"
lad<-"UCSC"
#list.files() will give us the names of the files and then we will iterate over them with sapply, open and calculate the mean
files <- list.files(path=work_dir, pattern="*.bed", full.names=TRUE, recursive=FALSE)

#Obtain a list with the mean length
mean_length<-sapply(files, function(x) {
  data <- as.data.frame(read.table(x, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  data$nt_intersecting<-(data$V3-data$V2)
  mean(data$nt_intersecting)
})

#Obtain a list with the median length
median_length<-sapply(files, function(x) {
  data <- as.data.frame(read.table(x, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  data$nt_intersecting<-(data$V3-data$V2)
  median(data$nt_intersecting)
})

#Obtain total number of bases
total_n_bases<-sapply(files, function(x) {
  data <- as.data.frame(read.table(x, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  data$nt_intersecting<-(data$V3-data$V2)
  sum(data$nt_intersecting)
})

#Iterate again to obtain then number of intersects. 
n_island<-sapply(files, function(x) {
  data <- as.data.frame(read.table(x, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  nrow(data)
})

#The data obtained from the loop is stored in lists, so we will convert into df 
mean_length_df<-as.data.frame(mean_length)
median_length_df<-as.data.frame(median_length)
total_n_bases_df<-as.data.frame(total_n_bases)
n_island_df<-as.data.frame(n_island)


#Then join both dataframes into 1
sum_table<-cbind(mean_length_df,median_length_df,total_n_bases_df, n_island_df)


#Transpose the df and change the name of the columns 
sum_table<-t(sum_table)
colnames(sum_table) <-c("H1.0","H1.2", "H1.4", "H1.5", "H1X")

sum_table<-t(sum_table)
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

sum_table<-round_df(sum_table, 0)
output_file<-paste(work_dir,cell_line,"_merged_islands.csv",sep="")
#Export the sum_table 
write.csv(sum_table,output_file)

