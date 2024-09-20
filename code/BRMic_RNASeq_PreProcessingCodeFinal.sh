#!/bin/bash

## Programs
trimgalore="/home/avm27/anaconda3/envs/RNA_Seq/bin/trim_galore"
cutadapt="/home/avm27/anaconda3/envs/RNA_Seq/bin/cutadapt"
STAR="/home/avm27/anaconda3/envs/RNA_Seq/bin/STAR"
SamTools="/home/avm27/anaconda3/envs/RNA_Seq/bin/samtools"
Fastqc="/home/avm27/anaconda3/envs/RNA_Seq/bin/fastqc"
Multiqc="/usr/bin/multiqc"
python="/home/avm27/anaconda3/envs/RNA_Seq/bin/python"
StringTie="/home/avm27/anaconda3/envs/RNA_Seq/bin/stringtie"
StringTieScript="/home/avm27/anaconda3/envs/RNA_Seq/bin/prepDE.py"
cores=7

## Experiment Details
experimentname="BRMicRNA"
experimentnamefull="Brain Region Microglia Experiment -- RNA Sequencing"


## Directories
maindir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq"
mergedfqdir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/MergedFQs"
trimmeddir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/TrimmedFQs"
stardir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/STAR_Output"
bamdir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/BAMs"
stringtiedir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/StringTie"
FPKMdir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/StringTie/FPKMs"
fastqcdir="/home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/FASTQC"

## Reference Files
mm10genome="/home/avm27/genome/mm10/ncbi_STAR"
mm10annoref="/home/avm27/genome/mm10/mm10.ncbiRefSeq.gtf"
mm10index="/home/avm27/genome/mm10/mm10.fa.fai"

echo "Final Code: BR Mic Experiment Pegasus Pre-Processing"

echo "Making Necessary Directories"

mkdir $mergedfqdir
mkdir $trimmeddir
mkdir $stardir
mkdir $bamdir
mkdir $stringtiedir
mkdir $FPKMdir
mkdir $fastqcdir

echo "Loading Conda Environment"

conda activate RNA_Seq

cd $maindir

## Concatenate R1 and R2 for each sample by lane

for sampleName in $(find $maindir -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
do 

    (echo "Merging R1 $sampleName"

    cat $maindir/${sampleName}_L00*_R1_001.fastq.gz > $mergedfqdir/${sampleName}_R1.fastq.gz
   
    echo "Merging R2 $sampleName"

    cat $maindir/${sampleName}_L00*_R2_001.fastq.gz > $mergedfqdir/${sampleName}_R2.fastq.gz) &

done

wait


echo "Trimming Samples"

for sampleName in $(find $mergedfqdir -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    $trimgalore --path_to_cutadapt $cutadapt --cores $cores --basename ${sampleName} --dont_gzip --paired $mergedfqdir/${sampleName}_R1.fastq.gz $mergedfqdir/${sampleName}_R2.fastq.gz --output_dir $trimmeddir &
    
done
wait


echo "Trimming Completed"

echo "Begin Aligning to the mm10 Genome"

for sampleName in $(find $mergedfqdir -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

   $STAR --runThreadN $cores \
    --runMode alignReads  \
    --genomeDir $mm10genome \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD XS \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNoverLmax 0.3 \
    --outFileNamePrefix $stardir/${sampleName} \
    --quantMode TranscriptomeSAM GeneCounts \
    --readFilesIn $trimmeddir/${sampleName}_val_1.fq $trimmeddir/${sampleName}_val_2.fq \
    --outTmpDir $stardir/${sampleName} &


done

wait

echo "Compressing Intermediate fq Files"

for sampleName in $(find $mergedfqdir -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    (gzip $trimmeddir/${sampleName}_val_1.fq
    gzip $trimmeddir/${sampleName}_val_2.fq) &

done

wait


echo "Obtaining unique mapped reads"

for sampleName in $(find $mergedfqdir -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    $SamTools view $stardir/${sampleName}Aligned.sortedByCoord.out.bam | grep -w 'NH:i:1' | $SamTools view -bt $mm10index > $bamdir/${sampleName}.bam &

done

wait


echo "Indexing Samples"

for sampleName in $(find $mergedfqdir -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    $SamTools index $bamdir/${sampleName}.bam &

done

wait

echo "Create GTF Files"

for sampleName in $(find $mergedfqdir -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    $StringTie -p $cores -G $mm10annoref -e -o $stringtiedir/${sampleName}.gtf -A $FPKMdir/${sampleName}_FPKM.txt -l ${sampleName} $bamdir/${sampleName}.bam &

done

wait

## Generate Sample Info List for Stringtie. Name the file "sample_lst.txt" and save to the working directory. Below is the code ran for this experiment.

echo "Generating sample info list for StringTie"

echo "MVTA6 /home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/StringTie/MVTA6.gtf" > $stringtiedir/sample_lst_$experimentname.txt
echo "MCTX4 /home/avm27/Documents/Raw_Sequencing_Data/BR_MIC_RNASeq/StringTie/MCTX4.gtf" >> $stringtiedir/sample_lst_$experimentname.txt


echo "Completed Same Info File"

## Generate Gene Count matrix for DESEQ2. Name the file "BRMicRNA_gene_count_matrix.csv" and save to the working directory.

echo "Generating Gene Count Matrix File"

$python $StringTieScript -i $stringtiedir/sample_lst_$experimentname.txt -g $stringtiedir/BRMicRNA_gene_count_matrix.csv

echo $experimentnamefull "RNA-Sequencing Pre-Processing Analysis Complete, Now completing FASTQC"

$Fastqc --outdir $fastqcdir -t 16 $maindir/*.fastq.gz 
$Fastqc --outdir $fastqcdir -t 16 $mergedfqdir/*.fastq.gz 
$Fastqc --outdir $fastqcdir -t 16 $trimmeddir/*

$Multiqc $maindir/ $mergedfqdir/ $trimmeddir/ $bamdir/ $stringtiedir/ $FPKMdir/


echo $experimentnamefull "QC Analysis Complete, Now Exiting"

