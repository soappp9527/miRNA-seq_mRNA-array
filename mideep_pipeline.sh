conda install -c bioconda mirdeep2
cd /data/2018_miRNA/

#trimming

#fastp
ls *_sequence.fastq.gz | parallel -j 2 fastp -i {} -o qc/{.}_out.fastq \
	-q 25 --length_limit 30 --adapter_fasta /data/adapters/TruSeq3-SE.fa

#trimmomatic
trimmomatic SE -threads 8 -phred33 ./1_sequence.fastq.gz ./qc/1_dd.fastq.gz \
	ILLUMINACLIP:/data/adapters/TruSeq3-SE.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
#bescause can't , need absolute path

#mapping
#mirdeep2

#build an index
bowtie-build -f /data/grch38/GRCh38.fasta grch38

#mapping
ls ./qc/*.fastq | parallel -j 2 mapper.pl {} -e -h -i -j -l 18 -m -p grch38 \
	-s ./mapping/{/}.fa -t ./mapping/{/}.arf -v -o 4

#extract known miRNA
extract_miRNAs.pl /data/mirbase/22.1/mature.fa hsa > mature_hsa.fa
extract_miRNAs.pl /data/mirbase/22.1/hairpin.fa hsa > hairpin_hsa.fa
extract_miRNAs.pl /data/mirbase/22.1/mature.fa mmu,chi > mature_other.fa
#no mention in offical document



#mirdeep2
miRDeep2.pl ./mapping/1_sequence.fastq_out.fastq.fa \
	genome.fa \
	./mapping/1_sequence.fastq_out.fastq.arf \
	mature_hsa.fa \
	mature_other.fa \
	hairpin_hsa.fa \
	-t hsa 2> report.log


#encounter error:
#Genome file /data/grch38/GRCh38.fasta has not allowed whitespaces in its first identifier
#need to remove everything after the space
#https://www.biostars.org/p/59181/
perl -plane 's/\s+.+$//' < /data/grch38/GRCh38.fasta > ./genome.fa

mkdir mideep
cd ./mideep

#
ls ../mapping/*.fa | parallel -j 1 miRDeep2.pl {.}.fa ../genome.fa {.}.arf ../mature_hsa.fa ../mature_other.fa ../hairpin_hsa.fa -t hsa 2> report.log


for x in ../mapping/*.fa
do
filename=$(basename "$x" .fa)
miRDeep2.pl ../mapping/$filename.fa ../genome.fa ../mapping/$filename.arf ../mature_hsa.fa ../mature_other.fa ../hairpin_hsa.fa -t hsa 2> report.log
done
