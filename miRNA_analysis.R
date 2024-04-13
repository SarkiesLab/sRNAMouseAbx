#miRNA_analysis

load("readCounts.Rdata")


library(DESeq2)


miRNAs<-read.table("miRNA_counts_all.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)

head(miRNAs)

#remove miRNAs with less than 5 rpm in any 

miRNA_normRPM<-miRNAs

for(i in 1:6){

miRNA_normRPM[,i]<-miRNA_normRPM[,i]*1e6/readCounts[i]



}

Filt<-which(apply(miRNA_normRPM,1,max)>5)
miRNA_filt<-miRNAs[Filt,]


sample_info<-cbind(1:6, c("C","C","C","T","T","T"))
colnames(sample_info)<-c("name","cond")


miRNA_deseq<-DESeqDataSetFromMatrix(miRNA_filt,sample_info, design=~cond)

miRNA_deseq<-DESeq(miRNA_deseq)

miRNA_deseqRes<-results(miRNA_deseq)

miRNA_deseqNorm<-counts(miRNA_deseq, normalized=T)

miRNA_deseqNorm<-cbind(miRNA_deseqNorm, apply(miRNA_deseqNorm[,1:3],1,mean),apply(miRNA_deseqNorm[,4:6],1,mean))

colnames(miRNA_deseqNorm)[7:8]<-c("mean_C","mean_T")

plot(log2(miRNA_deseqNorm[,7]),log2(miRNA_deseqNorm[,8]),col="gray", pch=16,xlab="log2_miRNAs_cont",ylab="log2_miRNAs_treat" )

sigMi<-which(miRNA_deseqRes$pvalue<0.05)

points(log2(miRNA_deseqNorm[sigMi,7]),log2(miRNA_deseqNorm[sigMi,8]),col="red", pch=16)
