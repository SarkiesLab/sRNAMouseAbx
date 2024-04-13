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



