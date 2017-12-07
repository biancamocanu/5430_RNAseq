#!/bin/sh

# Author: Bianca Mocanu
# MCB5430 RNAseq Final Project - shell pipeline
#===================================================================================================================================================

#This script contains all the Shell part of a RNAseq analysis pipeline starting from fastq files

#===================================================================================================================================================
# Requirements:

# 1.fastq files in inPATH
# 2.indexed genome or genome file to be indexed with "bowtie-build <infile> <outfile_handle>"
# 3.fastqc module (check for latest version: fastqc/0.11.5/, retrieved on Sep. 28, 2017)
# 4.fastx_tools (check for availability on BBC in /share/apps/ - add to $PATH if needed!
# 5.hisat2 (check for latest version: hisat2, retrieved on Dec. 6, 2017)
# 6.samtools (check for latest version: samtools/1.3.1/, retrieved on Oct. 19th, 2017)
# 7.bedtools (check for latest version: BedTools/2.26.0/, retrieved on Oct. 19th, 2017)
# 8. stringtie (check for latest version: stringtie/1.3.0, retrieved on Nov. 20th, 2017)
# 9. gffcompare

#===================================================================================================================================================
# Required modules load here:

module load bowtie2/2.3.1/
module load fastqc/0.11.5/
module load samtools/1.3.1/
module load BedTools/2.26.0/
module load stringtie/1.3.0/

#===================================================================================================================================================
# Global variables

inPATH="/tempdata3/MCB5430/midterm/midterm/fastq/" # uncomment this for the real (very large) files
# inPATH="/home/bim16102/midterm/" #used this on 1 mil reads files to test the script
hg19index="/tempdata3/MCB5430/genomes/hg19/bowtieIndex/hg19"
hg19chromInfo="/tempdata3/MCB5430/genomes/hg19/hg19_chromInfo.txt"
gencode="/tempdata3/MCB5430/annotations/hs/bed/hg19_gencode_ENSG_geneID.bed"
chr12="/tempdata3/MCB5430/genomes/hg19/fasta/chr12.fa"
hg19="/tempdata3/MCB5430/genomes/hg19/fasta/hg19.fa"
TSSbackground="/tempdata3/MCB5430/midterm/midterm/hg19_unique_TSSonly_bkgrnd.txt"
outPATH="/home/bim16102/final/processed_data/"
adapter="GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA"
summits_highconf=$(find /tempdata3/MCB5430/midterm/midterm/peaks/ -maxdepth 1 -type f)
fastqfiles=$(find ${inPATH} -maxdepth 1 -type f)
#===================================================================================================================================================

# Generating output directories

mkdir $outPATH
cd $outPATH
mkdir ./Qtrimmed_fastq
mkdir ./fastqc_reports
mkdir ./hisat_OUT
	mkdir ./hisat_OUT/sam
	mkdir ./hisat_OUT/bam
	mkdir ./hisat_OUT/bed
	mkdir ./hisat_OUT/bedgraphs
	
mkdir ./stringtie_OUT

for file in $fastqfiles
	do
		ext=`echo $(basename $file) | cut -d "." -f 2` # generated to see file type
		prefix=`echo $(basename $file) | cut -d "." -f 1`  #creates a prefix for each fastq file that is analyzed
		mkdir ./fastqc_reports/$prefix
		mkdir ./fastqc_reports/${prefix}_trimmed
		mkdir ./stringtie_OUT/$prefix
	
		if [ $ext=="fastq" ]
		then

			cd ./fastqc_reports/$prefix

			echo -e "Starting fastqc analysis on $(basename $file) ..."
			echo "Generating QC reports of unprocessed $(basename $file)"
			
			fastqc $file -o ./ 2>&1
			
			cd $outPATH/Qtrimmed_fastq/
			echo "Trimming low quality bases... (Q<30)"
			fastq_quality_trimmer -Q33 -t 30 -l 30 -i $file -o ./${prefix}_trimmed.fastq 2>&1

			echo "Generating QC reports of the quality trimmed $(basename $file)..."
			cd $outPATH/fastqc_reports/${prefix}_trimmed
			fastqc $outPATH/Qtrimmed_fastq/${prefix}_trimmed.fastq -o ./ 2>&1

			
			echo "Starting HISAT2 alignment to the hg19 gencode tran assembly"
			
		
			echo "Alignment to human genome (hg19) complete for $(basename $file)!"

			echo "Generating BAM file..."
			samtools view -S -b ${prefix}_chr12.sam > ${prefix}_chr12.bam  # sam to bam conversion

			samtools sort -l 9 -n  ${prefix}_chr12.bam -T ${prefix} -o ${prefix}_chr12.sorted.bam  #this sorts the bam file so that it occupies less space
			echo "BAM file generated!"

			echo "Generating BED file..."
			bedtools bamtobed -i ${prefix}_chr12.sorted.bam > ${prefix}_chr12.bed
			echo "BED file generated!"
			echo "Sorting BED file"
			sortBed -i ${prefix}_chr12.bed > ${prefix}_chr12_sorted.bed
			echo "Generating bedgraph file..."
			bedtools genomecov -ibam ${prefix}_chr12.sorted.bam -bg > ${prefix}.bedgraph #generates the bedgraph from bam directly
			echo "Bedgraph file generated!"

#==============================================================================================================================================================
# Comment this section if you want to keep all these files
#==============================================================================================================================================================
			echo "Cleaning up temporary files"
			rm ${prefix}_clipped.fastq # partially processed file
			rm ${prefix}_preprocessed.fastq # fully preprocessed fastq file - it occupies a lot of space
			rm ${prefix}.sam # large file as well, virtually useless once converted to bam
			rm ${prefix}_chr12.sam # subset of the file above, only for the assigned chromosome
			rm ${prefix}_chr12.bam # unsorted bam file
			rm ${prefix}_chr12.bed # unsorted bed file
		
		fi
#==============================================================================================================================================================

		echo "Preparing bedgraphs for Genome Browser..."
		
		if [ $prefix=="treatA_chip_rep1" ] || [ $prefix=="treatA_chip_rep2" ]
		then
			awk -v NAME="$prefix" 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=0,0,125"}
			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
		elif [ $prefix=="treatAB_chip_rep1" ] || [ $prefix=="treatAB_chip_rep2" ]
		then
			awk -v NAME=$prefix 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=125,0,125"}
			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
		elif [ $prefix=="Input" ]
		then
			awk -v NAME=$prefix 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=125,0,0"}
			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
		fi
		
		echo "Genome Browser bedgraphs generated!"
		
	done | tee -a ${outPATH}logfiles/log.txt