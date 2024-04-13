#readcounts

filelist<-list.files()
filelist<-filelist[grep("col_cut", filelist)]

readSum<-function(fileIn, minLength,maxLength){

fileTemp<-read.table(fileIn, sep="-", fill=T, stringsAsFactors=F)

fileTemp<-cbind(fileTemp[grep(">", fileTemp[,1]),],fileTemp[-grep(">", fileTemp[,1]),1])
print(head(fileTemp))
Filt<-which(nchar(fileTemp[,3])>=minLength&nchar(fileTemp[,3])<=maxLength)
print(length(Filt))


output<-sum(as.numeric(fileTemp[Filt,2]))
return(output)

}


readCounts<-c()
for(i in 1:6){

readCounts<-c(readCounts, readSum(filelist[i], 18,36))




}
names(readCounts)<-filelist
save(readCounts, file="readCounts.Rdata")