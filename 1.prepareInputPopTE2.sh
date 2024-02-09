#!/bin/bash

#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --job-name=bwa_waddington
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH -o logs/bwa_se2pe_%a.out
#SBATCH -e logs/bwa_se2pe_%a.err  

module load Miniconda3/4.7.10

source activate anopheles

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

#mkdir referenceGenome
#wget http://ftp.flybase.org/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-chromosome-r6.31.fasta.gz
#mv dmel-all-chromosome-r6.31.fasta.gz referenceGenome


# Create TE hierarchy file
echo -e "id\tfamily\torder" > $DIR/TE_annotation/TE_hierarchy.txt
for insertion in `grep ">" $DIR/TE_annotation/consensuses_curated_v4.fasta`
do
TEconsensus=$(echo $insertion | tr -d '>' )
order=$(grep -w "$TEconsensus" $DIR/TE_annotation/TE_library_family.csv | cut -f 6)
echo -e "$TEconsensus\t$TEconsensus\t$order" >> $DIR/TE_annotation/TE_hierarchy.txt
done

# mask reference genome
mkdir $DIR/referenceGenome/RepeatMasker

# heterochromatin from Rech et al 2022
echo -e "2L\t529999\t18870000\t2L\n2R\t5982494\t24972477\t2R\n3L\t749999\t19026900\t3L\n3R\t6754277\t31614278\t3R\nX\t1325966\t21338973\tX" > $DIR/referenceGenome/heterochromatin.bed

gunzip $DIR/referenceGenome/dmel-all-chromosome-r6.31.fasta.gz

#bedtools getfasta -fi $DIR/referenceGenome/dmel-all-chromosome-r6.31.fasta -bed $DIR/referenceGenome/heterochromatin.bed -nameOnly | seqtk seq -l60 > $DIR/referenceGenome/dmel-chr-eu-r6.31.fasta

echo -e "2L\n2R\n3L\n3R\nX" > $DIR/referenceGenome/list.txt 

seqtk subseq  $DIR/referenceGenome/dmel-all-chromosome-r6.31.fasta $DIR/referenceGenome/list.txt | seqtk seq -l60 > $DIR/referenceGenome/dmel-chr-r6.31.fasta

RepeatMasker -gccalc -s -cutoff 255 -no_is -nolow -norna -gff -u -pa 8 \
-lib  $DIR/TE_annotation/consensuses_curated_v4.fasta \
-dir $DIR/referenceGenome/RepeatMasker \
$DIR/referenceGenome/dmel-chr-r6.31.fasta

# concatenate masked genome with consensus TEs
mkdir -p $DIR/ANALYSES/PoPoolationTE2/input_PoPoolationTE2

cat $DIR/referenceGenome/RepeatMasker/dmel-chr-r6.31.fasta.masked $DIR/TE_annotation/consensuses_curated_v4.fasta > \
$DIR/ANALYSES/PoPoolationTE2/input_PoPoolationTE2/dmel.r6.31.TE.merged.fa

# Create index
bwa index $DIR/ANALYSES/PoPoolationTE2/input_PoPoolationTE2/dmel.r6.31.TE.merged.fa