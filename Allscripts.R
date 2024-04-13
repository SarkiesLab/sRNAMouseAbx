#readcounts
setwd("Col_fa")
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

#data_analysis

#piRNA analysis

load("readCounts.Rdata")


library(DESeq2)









#load piRNA data

piRNAs<-read.table("piRNA_count_clusters_all_min24_max36_fNT.txt", sep="\t", stringsAsFactors=F,header=T, row.names=1)




#filter for expressed piRNA clusters
#use crude normalization (rpm) to identify those which have expression above 0.1 rpm in any sample
piRNAs_normRPM<-piRNAs
for(i in 1:6){piRNAs_normRPM[,i]<-piRNAs_normRPM[,i]*1e6/readCounts[i]}
piRNAs_filt<-piRNAs[which(apply(piRNAs_normRPM,1,max)>=0.1),]
piRNAs_normRPM_filt<-piRNAs_normRPM[Filt,]
piRNA_sums<-apply(piRNAs_normRPM_file,2,sum)

sample_info<-cbind(1:6, c("C","C","C","T","T","T"))
colnames(sample_info)<-c("name","cond")



piRNA_deseq<-DESeqDataSetFromMatrix(piRNAs_filt,sample_info, design=~cond)

piRNA_deseq<-DESeq(piRNA_deseq)

piRNA_deseqRes<-results(piRNA_deseq)
#there are some with identical counts, duplicates should be removed as they are overlapping annotations
Dups<-which(duplicated(piRNA_deseqRes[,1])&duplicated(piRNA_deseqRes[,2]))
piRNA_deseqRes<-piRNA_deseqRes[-Dups,]
piRNA_normCounts<-counts(piRNA_deseq, normalized=T)
piRNA_normCounts<-piRNA_normCounts[-Dups,]

piRNA_normCounts<-cbind(piRNA_normCounts, apply(piRNA_normCounts[,1:3],1,mean), apply(piRNA_normCounts[,4:6],1,mean))
colnames(piRNA_normCounts)[7:8]<-c("mean_C","mean_T")

plot(log2(piRNA_normCounts[,7]),log2(piRNA_normCounts[,8]),pch=16,col="gray", ylab="log2_piRNA_clusters_Cont",xlab="log2_piRNA_clusters_Test")

sigPi<-which(piRNA_deseqRes$pvalue<0.05)
points(log2(piRNA_normCounts[sigPi,7]),log2(piRNA_normCounts[sigPi,8]),col="red", pch=16)


legend("topleft", "p<0.05",pch=16, col="red")

dev.copy(pdf, "piRNA_cluster_comparison.pdf")
dev.off()
piRNA_sums<-apply(piRNAs_normRPM_filt,2,sum)
save(piRNA_sums, file="piRNA_sums.Rdata")
setwd("../")
#tRNA analysis
#tRNA analysis
#individual small RNA sequences from control
#these are mapped to tRNAs using bedtools intersectBed 
#this file acts as the input for a differential expression analysis for 5'tRNAs and 3'tRNAs


#in R
setwd("tRNA_analysis")
library(DESeq2)

tRNA_annotation<-read.table("tRNAs_seq.bed", sep="\t", stringsAsFactors=F)

#now count these individual sequences in each file separately

tRNA_seqs<-tRNA_annotation[,9]

tRNA_files<-list.files()

tRNA_counts<-matrix(0, nrow=length(tRNA_seqs), ncol=6)

for(i in 1:6){

#read in file
temp<-read.table(tRNA_files[i], sep="\t", stringsAsFactors=F)

#extract readcounts

count_temp<-read.table(text=temp[,4], sep="-", stringsAsFactors=F)
count_temp<-as.numeric(count_temp[,2])

#match by sequence
seqMatch<-match(tRNA_seqs, temp[,9])

#populate counts
tRNA_counts[,i]<-count_temp[seqMatch]
tRNA_counts[-which(tRNA_counts[,i]>0),i]<-0
temp<-c()
count_temp<-c()
seqMatch<-c()
}


colnames(tRNA_counts)<-tRNA_files

row.names(tRNA_counts)<-paste(tRNA_annotation[,9], tRNA_annotation[,13], sep="|")

#separate into 5' and 3' halves

sense_tRNAs_pos<-which(tRNA_annotation[,6]=="+"&tRNA_annotation[,15]=="+")

sense_tRNAs_neg<-which(tRNA_annotation[,6]=="-"&tRNA_annotation[,15]=="-")

Start_dist<-rep(-100, length=nrow(tRNA_annotation))
#populate this with distance between start of sRNA seq and start of tRNA 
Start_dist[sense_tRNAs_pos]<-tRNA_annotation[sense_tRNAs_pos,2]-tRNA_annotation[sense_tRNAs_pos,16]
Start_dist[sense_tRNAs_neg]<-tRNA_annotation[sense_tRNAs_neg,3]-tRNA_annotation[sense_tRNAs_neg,17]
#select sequences that start within 1 bp of the tRNA start and are short <35bp

tRNA5s<-which(abs(Start_dist)<1&tRNA_annotation[,7]>=30&tRNA_annotation[,7]<36)

tRNA5_counts<-tRNA_counts[tRNA5s,]


#do the same for the end of the sequence to map tRNA3s
End_dist<-rep(-100, length=nrow(tRNA_annotation))
End_dist[sense_tRNAs_pos]<-tRNA_annotation[sense_tRNAs_pos,3]-tRNA_annotation[sense_tRNAs_pos,17]
End_dist[sense_tRNAs_neg]<-tRNA_annotation[sense_tRNAs_neg,2]-tRNA_annotation[sense_tRNAs_neg,16]


tRNA3s<-which(abs(End_dist)<1&tRNA_annotation[,7]<36&tRNA_annotation[,7]>=30)


tRNA3_counts<-tRNA_counts[tRNA3s,]

#now we have two separate counts tables.

#ready for deseq
sample_info<-data.frame(1:6, c("C","C","C","T","T","T"))
#where C is control (untreated) and T is treated

colnames(sample_info)<-c("name","cond")

tRNA5_deseq<-DESeqDataSetFromMatrix(tRNA5_counts, sample_info, design=~cond)
tRNA5_deseq<-DESeq(tRNA5_deseq)
tRNA5_deseqRes<-results(tRNA5_deseq)

normtRNA5<-counts(tRNA5_deseq, normalized=T)
normtRNA5<-cbind(normtRNA5, apply(normtRNA5[,1:3],1,mean), apply(normtRNA5[,4:6],1,mean))

colnames(normtRNA5)[7:8]<-c("mean_C","mean_T")
plot(log2(normtRNA5[,7]), log2(normtRNA5[,8]),pch=16,col="gray",xlab="log2_5tRNA_cont",ylab="log2_5tRNA_treat")
sigs<-which(tRNA5_deseqRes$pvalue<0.05)
points(log2(normtRNA5[sigs,7]),log2(normtRNA5[sigs,8]), pch=16,col="red")
legend("topleft","p<0.05", pch=16,col="red" )
text(18,20, "tRNA_Gly_GCC")
dev.copy(pdf, "tRNA_5_deseq.pdf")
dev.off()

tRNA3_deseq<-DESeqDataSetFromMatrix(tRNA3_counts, sample_info, design=~cond)
tRNA3_deseq<-DESeq(tRNA3_deseq)
tRNA3_deseqRes<-results(tRNA3_deseq)

normtRNA3<-counts(tRNA3_deseq, normalized=T)
normtRNA3<-cbind(normtRNA3, apply(normtRNA3[,1:3],1,mean),apply(normtRNA3[,4:6],1,mean))

colnames(normtRNA3)[7:8]<-c("mean_C","mean_T")

plot(log2(normtRNA3[,7]),log2(normtRNA3[,8]),pch=16, col="gray", xlab="log2_3tRNA_cont", ylab="log2_3tRNA_treat")
sigs3<-which(tRNA3_deseqRes$pvalue<0.05)
points(log2(normtRNA3[sigs3,7]),log2(normtRNA3[sigs3,8]),col="red",pch=16)
legend("topleft", "p<0.05",pch=16,col="red")
text(3,6, "tRNA-Asp-GTC")
dev.copy(pdf, "tRNA_3_deseq.pdf")
dev.off()
tRNA5_sums<-apply(tRNA5_counts,2,sum)
tRNA3_sums<-apply(tRNA3_counts,2,sum)

save(tRNA5_sums, file="tRNA5_sums.Rdata")
save(tRNA3_sums, file="tRNA3_sums.Rdata")

#specific look at gly-tcc

boxplot(normtRNA5[3,1:3],normtRNA5[3,4:6], names=c("Control","Treat"), ylab="normalized_reads_deseq",col="white")
points(c(1,1,1,2,2,2),normtRNA5[3,1:6],col=c("blue","blue","blue","purple","purple","purple"), pch=16)
setwd("../")
#miRNA_analysis

#miRNA_analysis
setwd("miRNA_analysis")
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

