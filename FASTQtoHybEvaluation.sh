#!/bin/bash

############################################################################
# User Inputs Below
############################################################################

#OPTIONAL: Change ONLY if not executing in same directory as fastq files. 
origin=$(pwd)

referencegenome = 'absolute path to bowtie indexes'

referencefasta = 'absolute path to reference genome fa/fasta'

PanelName = 'name of panel'

baitintervalfile = 'absolute path to bait interval list'

targetintervalfile = 'absolute path to target interval list'


#end of user inputs
########################################################################


#check if downsampling is required, and to what number of CLUSTERS
downsample=''
while [[ $downsample != 'Y' && $downsample != 'N' ]]
do
  echo 'Downsample? (Y/N)'
  read downsample
  if [ $downsample == 'Y' ]; then
    #get number of clusters
    echo "How many CLUSTERS? remember for paired end sequencing 1 cluster = 2 reads"
    read clusters
  elif [ $downsample == 'N' ]; then
    #if no downsampling, process all the reads in each file
    clusters='allreads'
  else
    #didn't input Y or N, start loop over
    echo 'Not a valid input'
  fi
done


#Define Merge Function
Merge () {
  #'''Merges multiple fastq files for same sample if multiple lanes'''
  #'''$1 = Read (R1 or R2)'''
  pwd_dir=$(pwd)
  pwd_fastq_R="${pwd_dir}/*L002_${1}_001.fastq.gz"
  for sample in `ls ${pwd_fastq_R}`
  do
    dir=${pwd_dir}
    base=$(basename $sample "L002_${1}_001.fastq.gz")
    echo "Merging $base $1"
    cat ${dir}/${base}L001_${1}_001.fastq.gz ${dir}/${base}L002_${1}_001.fastq.gz ${dir}/${base}L003_${1}_001.fastq.gz ${dir}/${base}L004_${1}_001.fastq.gz > ${dir}/M_${base}${1}_001.fastq.gz
    rm ${dir}/${base}L00*_${1}_001.fastq.gz
    #QC fastq file
    if [ ${1} == "R1" ]; then
      fastqc ${dir}/M_${base}${1}_001.fastq.gz
    fi
  done
}


#Do trimmed fastq files already exist? if yes, skip merging
# and trimming, if not, merge and trim and then continue on
# allows for faster reanlysis
DIR="${origin}/Trim"
# check if trim directory already exists, and if it does, does it have files in it?
if [ -d "$DIR" ]; then
	if [ "$(ls -A $DIR)" ]; then
    cd Trim
    fi
else 
    #merge if it hasn't been done yet
    Merge R1
    Merge R2

    #Trim adapters
    pwd_dir=$(pwd)
    pwd_fastq="${pwd_dir}/*R1_001.fastq.gz"
    mkdir Trim
    for sample in `ls ${pwd_fastq}`
    do
      dir=${pwd_dir}
      base=$(basename $sample "R1_001.fastq.gz")
      echo "File 1: ${dir}/${base}R1_001"
      echo "File 2: ${dir}/${base}R2_001"

      time trimmomatic PE ${dir}/${base}R1_001.fastq.gz ${dir}/${base}R2_001.fastq.gz \
      ${dir}/Trim/${base}Trimmed_R1_001.fq.gz ${dir}/${base}UnPaired_R1_001.fq.gz \
      ${dir}/Trim/${base}Trimmed_R2_001.fq.gz ${dir}/${base}UnPaired_R2_001.fq.gz \
      ILLUMINACLIP:TruSeq3-PE-2.fa:1:10:5:9:true MINLEN:50
    done

    #remove unpaired reads
    rm *UnPaired_R*_001.fq.gz

    cd Trim

    #just R2 this time, QC trimmed fastq files
    fastqc *R2_001_fq.gz
fi


#downsampling if relevant
if [ $clusters != 'allreads' ]; then
  #downsample
  pwd_dir=$(pwd)
  pwd_fastq_trimmed="${pwd_dir}/*Trimmed_R1_001.fq.gz"
  mkdir ${origin}/Subsample_${clusters}
  for sample in ${pwd_fastq_trimmed}
  do
    dir=${pwd_dir}
    base=$(basename $sample "Trimmed_R1_001.fq.gz")
    echo "Downsampling ${sample} to ${clusters} clusters"
    seqtk sample -s100 ${dir}/${base}Trimmed_R1_001.fq.gz ${clusters} | gzip > ${origin}/Subsample_${clusters}/Sub_${base}Trimmed_R1_001.fq.gz
    seqtk sample -s100 ${dir}/${base}Trimmed_R2_001.fq.gz ${clusters} | gzip > ${origin}/Subsample_${clusters}/Sub_${base}Trimmed_R2_001.fq.gz
  done
  echo "Downsampling done"
  cd ${origin}/Subsample_${clusters}
fi


#Alignment with Bowtie2
pwd_dir=$(pwd)
#Align with bowtie2
pwd_fastq_trimmed="${pwd_dir}/*Trimmed_R1_001.fq.gz"
mkdir ${origin}/Alignment_${clusters}
#Align with bowtie2
for sample in ${pwd_fastq_trimmed}
do
  dir=${pwd_dir}
  base=$(basename $sample "Trimmed_R1_001.fq.gz")
  echo "File 1: ${base}"
  echo "File 2: ${base}"
  echo "Outputfile: ${base}bowtie2"
  #do alignment with bowtie2 and create log of alignment efficiency
  time bowtie2 -x $referencegenome \
  -1 ${dir}/${base}Trimmed_R1_001.fq.gz -2 ${dir}/${base}Trimmed_R2_001.fq.gz \
  -p 4 -S ${origin}/Alignment_${clusters}/${base}Bowtie2.sam 2> ${origin}/Alignment_${clusters}/${base}Bowtie2.log
done

echo "Alignment done"

cd ${origin}/Alignment_${clusters}

#convert sam files to sorted bam files
for file in *.sam
do
  echo "Sorting file: $file"; time samtools sort -O bam -o ${file%_*}_sorted.bam $file
  echo "done"
done

#mark duplicates and output marked file
for file in *_sorted.bam
do
  echo "Marking duplicate: $file"
  time picard MarkDuplicates \
  I=$file \
  O=${file%_sorted.bam}_Marked.bam \
  M=${file%_sorted.bam}_dup_metrics_Marked.txt \
  REMOVE_DUPLICATES=FALSE
done

#remove duplicates and output removed file
for file in *_sorted.bam
do
  echo "Removing duplicate: $file"
  time picard MarkDuplicates \
  I=$file \
  O=${file%_sorted.bam}_Removed.bam \
  M=${file%_sorted.bam}_dup_metrics_Removed.txt \
  REMOVE_DUPLICATES=TRUE
done

rm *_sorted.bam

#index bam files
for file in *.bam; do echo "Indexing file: $file"; time samtools index $file; done


#define HSMetrics function
HS_Metrics () {
  #'''performs HS metrics for the indicated target file and compiles results
  # 1- PanelName
  # 2- Probe file path
  # 3- Target file path'''

  echo "Collecting ${1} HS_Metrics"

  #Make directories
  mkdir HSMetrics_${1}_${clusters}
  mkdir IntervalCoverage_${1}_${clusters}

  #Perform Picard HSMetrics
  for file in *.bam
  do
    base=$(basename $file ".bam")
    echo "collecting ${1} HS Metrics for $file"
    picard CollectHsMetrics \
    -I $file \
    -O ./HSMetrics_${1}_${clusters}/${base}_${1}_metrics.txt \
    -R $referencefasta \
    --BAIT_INTERVALS ${2} \
    --TARGET_INTERVALS ${3} \
    --PER_TARGET_COVERAGE ./IntervalCoverage_${1}_${clusters}/${file%_L00*}.intervals.txt 
  done

  cd HSMetrics_${1}_${clusters}

  #Extract HS Metrics
  today=$(date '+%Y_%m_%d')
  operator=$USER
  echo "generated on" $today "by" $operator >> HSMetrics.tsv
  echo "Probe File" $2 >> HSMetrics.tsv
  echo "Target File" $3 >> HSMetrics.tsv
  for files in *.txt
    do
    cat $files | awk 'NR==7, OFS="\t" { print "FileName", $0}' >> HSMetrics.tsv
    break
    done
    for files in *.txt
    do
    cat $files | awk -v temp="${files%_metrics.txt}" 'NR==8, OFS="\t" {print temp, $0}' >> HSMetrics.tsv
  done

  sed 's/\t/,/g' HSMetrics.tsv > $(date +"%Y%m%d")_HSMetrics_${1}_${clusters}.csv
  rm HSMetrics.tsv

  echo "${1} HS Metrics Done"

  cd ..
}


#Run HS_Metrics, uses provided panel info
HS_Metrics $panelname $baitintervalfile $targetintervalfile


#collect insert size Metrics
mkdir InsertSize_${clusters}
#collect insert size metrics
for sample in *Marked.bam
do
  picard CollectInsertSizeMetrics \
  I=$sample \
  O=./InsertSize_${clusters}/${sample}_InsertSize_Marked.txt \
  H=./InsertSize_${clusters}/${sample}_InsertSizeHistogram_Marked.pdf
done


cd ${origin}/Alignment_${clusters}

rm *Bowtie2.sam

#Final QC
cd ${origin}

multiqc . -n MultiQC_${today}.html --fn_as_s_name -f

echo 'Analysis Finished'
