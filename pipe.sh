#!/bin/sh

# Author: Bianca Mocanu
# MCB5430 RNAseq Final Project - shell pipeline
#===================================================================================================================================================

#This script contains all the Shell part of a RNAseq analysis pipeline starting from fastq files

#===================================================================================================================================================
# Requirements:

# 1.fastq files in inPATH
# 2.annotated genome for mapping with HISAT2 - not doable on your laptop, but you can run scripts that require tons of memory (160 GB ram for >2 hrs) - or download the already annotated genome from the HISAT2 website, if available
# 3.fastqc module (check for latest version: fastqc/0.11.5/, retrieved on Sep. 28, 2017)
# 4.fastx_tools (check for availability on BBC in /share/apps/ - add to $PATH if needed!
# 5.hisat2 (check for latest version: hisat2, retrieved on Dec. 6, 2017)
# 6.samtools (check for latest version: samtools/1.3.1/, retrieved on Oct. 19th, 2017)
# 7.bedtools (check for latest version: BedTools/2.26.0/, retrieved on Oct. 19th, 2017)
# 8.stringtie (check for latest version: stringtie/1.3.0, retrieved on Nov. 20th, 2017)
# 9.prepDE python script to create the files to be loaded in RStudio later (/tempdata3/MCB5430/scripts/prepDE.py)

#===================================================================================================================================================
# Required modules load here:

module load fastqc/0.11.5/
module load BedTools/2.26.0/
module load stringtie/1.3.0/
module load hisat2

#===================================================================================================================================================
# Global variables

inPATH="/archive/MCB5430/final_data/fastq/"
hg19index="/tempdata3/MCB5430/genomes/hg19/hisat_gencode_tran/hg19_gencode"
hg19chromInfo="/tempdata3/MCB5430/genomes/hg19/hg19_chromInfo.txt"
gencode="/tempdata3/MCB5430/annotations/hs/bed/hg19_gencode_ENSG_geneID.bed"
gencode_GTF="/tempdata3/MCB5430/annotations/hs/GTF/gencode.v19.annotation_pc.gtf"
hg19="/tempdata3/MCB5430/genomes/hg19/fasta/hg19.fa"
outPATH="/tempdata3/MCB5430/Bianca_Final/processed_data/"  #using this folder since there's a size limit for the individual accounts' storage
fastqfiles=$(find ${inPATH} -maxdepth 1 -type f)
prepDE="/tempdata3/MCB5430/scripts/prepDE.py"
#===================================================================================================================================================
# Generating output directories

mkdir /tempdata3/MCB5430/Bianca_Final/
mkdir $outPATH
cd $outPATH
mkdir ./Qtrimmed_fastq
mkdir ./fastqc_reports
mkdir ./hisat_OUT
	mkdir ./hisat_OUT/sam
	mkdir ./hisat_OUT/bam
	mkdir ./hisat_OUT/bedgraphs
mkdir ./logfiles
mkdir ./stringtie_OUT

#====================================================================================================================================================

for file in $fastqfiles
	do
		cd $outPATH
		ext=`echo $(basename $file) | cut -d "." -f 2` # generated to see file type
		prefix=`echo $(basename $file) | cut -d "." -f 1`  #creates a prefix for each fastq file that is analyzed
		mkdir ./fastqc_reports/$prefix
		mkdir ./fastqc_reports/${prefix}_trimmed
		mkdir ./stringtie_OUT/$prefix
	
		if [ $ext=="fastq" ]
		then

#====================================================================================================================================================
# The steps for fastQC of the initial files, quality trimming (for Q<30 and the remaining fragment length <30), and fastQC on the trimmed files go below:		
	
			cd ./fastqc_reports/$prefix 

			echo -e "Starting fastqc analysis on $(basename $file) ..."
			echo "Generating QC reports of unprocessed $(basename $file)"
			
			fastqc $file -o ./ 2>&1 
			touch ${prefix}_fastqcgoeshere.txt
			
			cd $outPATH/Qtrimmed_fastq/
			echo "Trimming low quality bases... (Q<30)"
			fastq_quality_trimmer -Q33 -t 30 -l 30 -i $file -o ./${prefix}_trimmed.fastq 2>&1

			echo "Generating QC reports of the quality trimmed $(basename $file)..."
			cd $outPATH/fastqc_reports/${prefix}_trimmed
			fastqc $outPATH/Qtrimmed_fastq/${prefix}_trimmed.fastq -o ./ 2>&1
			touch ${prefix}_fastqc_trimmed.txt
			
# HISAT 2 alignment of the quality trimmed files - use 4 processors and the given hg19 hisat index

			echo "Starting HISAT2 alignment to the hg19 gencode tran assembly"
			cd $outPATH/hisat_OUT/sam
			hisat2 -p4 --dta -x $hg19index -U $outPATH/Qtrimmed_fastq/${prefix}_trimmed.fastq -S ./${prefix}.sam  2>&1
			
			echo "Alignment to human genome (hg19) complete for $(basename $file)!"

			echo "Generating BAM file..."
			samtools view -S -b ${prefix}.sam > ../bam/${prefix}.bam  # sam to bam conversion  2>&1
			cd ../bam
			samtools sort -l 9 -n ${prefix}.bam -T ${prefix} -o ${prefix}_sorted.bam -g $hg19_chromInfo 2>&1
			
			rm ${prefix}.bam 
			
			echo "BAM file generated!"

# The genomecov line below takes FOREVER to run on the cluster		
			echo "Generating bedgraph files for browser data visualization..."
			bedtools genomecov -bga -split -ibam ./${prefix}_sorted.bam -g $hg19chromInfo > ../bedgraphs/${prefix}_IGV.bedgraph 2>&1
			
			echo "Preparing bedgraphs for Genome Browser..."
		
# I know it's hard coded but I haven't found a single shell script so far that doesn't have "filename rules" in the description so that they do it better.

			if [ $prefix=="E2_rep1" ] || [ $prefix=="E2_rep2" ]
			then
				awk -v NAME="$prefix" 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
				print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=0,0,125"}
				{print $0}' ${prefix}_IGV.bedgraph > ${prefix}_UCSC.bedgraph
			else
				awk -v NAME=$prefix 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
				print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=125,0,125"}
				{print $0}' ${prefix}_IGV.bedgraph > ${prefix}_UCSC.bedgraph
			fi
		
			echo "Genome Browser bedgraphs generated!"
			
			echo "Bedgraph file generated!"
			
# This is the stringtie part of the script - stringtie by default takes a bam file as input. The final asks us to use .sam as input, but sorting the sam file outputs a bam file regardless. I'm therefore using the sorted bam output in stringtie because I can't be bothered to find a fix, and it works.
			
			samtools sort -@ 4 ${prefix}.sam ${prefix}_sorted
			stringtie -e -B -p 4 -G $gencode_GTF -o ../../stringtie_OUT/${prefix}_est.gtf -A ../../stringtie_OUT/${prefix}_abun.tab ./${prefix}.sam

# This creates the sample_list file for the prepDE.py script. Important to note that if you run this script twice and you have 4 samples, this sample_list file will have 8 files in it. I might do a cleanup at the end that delets this file because it's relatively useless after prepDE anyway.
 
			cd ../../stringtie_OUT
			echo "$prefix	${prefix}_est.gtf" >> ./sample_list.txt
						
		fi
# Command for the prepDE script using the sample_list generated above		
		cd $outPATH/stringtie_OUT/
		python $prepDE -i sample_list.txt
			
	done | tee -a ${outPATH}/logfiles/log.txt
