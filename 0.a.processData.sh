#!/bin/bash
#SBATCH --mem=16G
#SBATCH --job-name=data_processing
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH -o logs/process_data.out 
#SBATCH -e logs/process_data.err 

DIR=/homes/users/mcoronado/scratch/Waddington/

for folder in `ls $DIR/DATA`
do
while read fastq_md5
do
md5=$(echo "$fastq_md5" | cut -f1 -d' ')
file=$(echo "$fastq_md5" | cut -f3 -d' ')
md5download=$(md5sum $DIR/DATA/$folder/$file | cut -f1 -d' ')
if [[ "$md5" == "$md5download" ]]
then
echo "$folder $file $md5 $md5download ok"
else
echo "$folder $file $md5 $md5download problem"
fi
done < $DIR/DATA/$folder/MD5.txt
done

for folder in `ls $DIR/DATA`
do
nFiles=$(ls $DIR/DATA/$folder | grep gz$ | wc -l)
files=$(ls $DIR/DATA/$folder | grep gz$ )

> $DIR/samples.tab
while read file
do
echo -e "$folder\t$file\t$nFiles"
done <<< "$files"
done >> $DIR/samples.tab

for folder in `ls $DIR/DATA`
do
mkdir -p $DIR/DATA_processed/$folder
rep=$(awk -v folder="$folder" ' $1 == folder ' samples.tab | cut -f 3 | sort -u )
if [[ "$rep" -eq 4 ]]
then
#echo $folder $rep
cat $DIR/DATA/$folder/*_1.fq.gz > $DIR/DATA_processed/$folder/${folder}_1.fq.gz
cat $DIR/DATA/$folder/*_2.fq.gz > $DIR/DATA_processed/$folder/${folder}_2.fq.gz
elif [[ "$rep" -eq 2 ]]
then
#echo $folder $rep
cp $DIR/DATA/$folder/*_1.fq.gz $DIR/DATA_processed/$folder/${folder}_1.fq.gz
cp $DIR/DATA/$folder/*_2.fq.gz $DIR/DATA_processed/$folder/${folder}_2.fq.gz
fi
done

for folder in `ls $DIR/DATA_processed`
do
echo "$folder" >> samples_processed.lst
done
