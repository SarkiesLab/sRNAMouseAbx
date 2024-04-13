export PATH=$PATH:/Users/psarkies/bowtie/
export PATH=$PATH:/Users/psarkies/samtools-1.9/
export PATH=$PATH:/Users/psarkies/bedtools2/bin/

#bowtie-build /Users/psarkies/bowtie/genomes/Mus_musculus.GRCm38.dna_sm.toplevel.fa Mus_GRCm38

#for i in `ls *.fa`
#do
#bowtie -f --sam -v 0 -k 1 Mus_GRCm38 "$i" >Aln_"$i".sam
#samtools view -bS Aln_"$i".sam >Aln_"$i".bam
#bamToBed -i Aln_"$i".bam >Aln_"$i".bed 
#rm Aln_"$i".sam
#done

#extract sequence information
#perl ../../../extract_seqinfo_bed.pl 

#extract piRNAs

#for i in `ls *.bed`
#do
#intersectBed -a "$i" -b ../piRNA_clusters_mm10.bed -wao |grep -v 0 >piRNA_"$i"
#done
#extract miRNAs
for i in `ls *.bed`
do
intersectBed -a "$i" -b /Users/psarkies/Utils/mmu_mirs_2.gff3 -wao | grep -v 0 >miRNA_"$i"
done
