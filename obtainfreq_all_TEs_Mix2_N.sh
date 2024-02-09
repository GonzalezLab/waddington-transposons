#!/bin/bash

#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --job-name=obtain_freq
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH -o logs/obtain_freq_%a.out
#SBATCH -e logs/obtain_freq_%a.err  
#SBATCH --exclude=mr-03-[01-26],mr-03-28,mr-04-00

module load Miniconda3/4.7.10

source activate anopheles

module load R/3.6.0-foss-2018b

# Run

# sbatch --array=1-1048%50 obtainfreq_all_TEs_Mix2_N.sh

DIR="/homes/users/mcoronado/scratch/Waddington"

cd $DIR

#types="A N"
#classes=$(cut -f1  $DIR/samples_processed.lst | cut -f1 -d'_' | sort -u)


#data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat ${DIR}/ANALYSES/post_analysis/data_comb.tab ))
class=Mix2
type=N

cd $DIR/ANALYSES/post_analysis/
# > $DIR/ANALYSES/post_analysis/frequency_insertions_all_${class}_$type.tab
# > $DIR/ANALYSES/post_analysis/fisher_PopTE2_all_${class}_$type.tab
# > $DIR/ANALYSES/post_analysis/fisher_TEMP2_all_${class}_$type.tab
# > $DIR/ANALYSES/post_analysis/FC_PopTE2_all_${class}_$type.tab
# > $DIR/ANALYSES/post_analysis/FC_TEMP2_all_${class}_$type.tab


#bedtools window -a ${DIR}/ANALYSES/common_insertions/reliable/${class}_${type}/${class}_${type}.teinsertion.eu.popte2.freq.temp2.reliable.bed -b ${DIR}/referenceGenome/dmel-all-r6.31.gtf > ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_all.intersect_1kb.bed

if [[ $class == "Mix2" ]]; then
class2="Mix"
elif [[ $class == "Mix" ]]; then
class2="Mix"
elif [[ $class == "D907" ]]; then
class2="D907"
fi

#nTes=$(cat  ${DIR}/ANALYSES/common_insertions/reliable/${class}_${type}/${class}_${type}.teinsertion.eu.popte2.freq.temp2.reliable.bed | wc -l )
TEinfo=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat  ${DIR}/ANALYSES/common_insertions/reliable/${class}_${type}/${class}_${type}.teinsertion.eu.popte2.freq.temp2.reliable.bed ))

#TEinfo=$(cat  ${DIR}/ANALYSES/common_insertions/reliable/${class}_${type}/${class}_${type}.teinsertion.eu.popte2.freq.temp2.reliable.bed )
#nGenes=$(cut -f17 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_all.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u | wc -l)
#genes=$(cut -f17 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_all.intersect_1kb.bed | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u | wc -l)
#wingUR=$(comm -12 <(cut -f17 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_all.intersect_1kb.bed  | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_wing_${class2}_${type}.lst | sort -u ) | wc -l )

#echo "$class $type - wing UR"
#genes=$(comm -12 <(cut -f17 ${DIR}/ANALYSES/post_analysis/genes/${class}_${type}_all.intersect_1kb.bed  | cut -f1 -d';' | cut -f 2 -d' ' | tr -d '"' | sort -u) <(cat $DIR/list_genes_DE/genes_UR_wing_${class2}_${type}.lst | sort -u ))

#while read TEinfo
#do
chr=$(echo "$TEinfo" | cut -f 1)
start=$(echo "$TEinfo" | cut -f 2)
end=$(echo "$TEinfo" | cut -f 3)
family=$(echo "$TEinfo" | cut -f 6)
TE=$(echo "$TEinfo" | cut -f 4)

# freqs=$(awk -v start="$start" ' $2 == start ' $DIR/ANALYSES/common_insertions/all.teinsertion.eu.popte2.freq.temp2.reliable.bed | awk -v TE="$TE" ' $4 == TE' | cut -f 10)

# freqC=$(echo $freqs | tr ';' '\n' | grep "${class}_C" | cut -f2 -d'=' | cut  -f1 -d';')
# freqA=$(echo $freqs | tr ';' '\n' | grep "${class}_A" | cut -f2 -d'=' | cut  -f1 -d';')
# freqN=$(echo $freqs | tr ';' '\n' | grep "${class}_N" | cut -f2 -d'=' | cut  -f1 -d';')
# freqP0=$(echo $freqs | tr ';' '\n' | grep "${class}_P0" | cut -f2 -d'=' | cut  -f1 -d';')

echo -e "${class}_$type\t$TE\t$start\t$end"

#echo -e "${class}_$type\twing_UR\t$TE\tPopTE2\t$freqP0\t$freqC\t$freqA\t$freqN" >> $DIR/ANALYSES/post_analysis/frequency_insertions_all_${class}_$type.tab

echo -e "$chr\t$start\t$end\t$TE\t$family" > $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed

Rscript family.R $DIR/ANALYSES/PoPoolationTE2/output/${class}_C/${class}_C.teinsertions
freqC=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_C/${class}_C.teinsertions.family) | grep -E "($family.*){2}"  |  cut -f 13)

Rscript family.R $DIR/ANALYSES/PoPoolationTE2/output/${class}_A/${class}_A.teinsertions
freqA=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_A/${class}_A.teinsertions.family) | grep -E "($family.*){2}"  |  cut -f 13)

Rscript family.R $DIR/ANALYSES/PoPoolationTE2/output/${class}_N/${class}_N.teinsertions
freqN=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_N/${class}_N.teinsertions.family) | grep -E "($family.*){2}"  |  cut -f 13)

Rscript family.R $DIR/ANALYSES/PoPoolationTE2/output/${class}_P0/${class}_P0.teinsertions
freqP0=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_P0/${class}_P0.teinsertions.family) | grep -E "($family.*){2}"  |  cut -f 13)

nFreqC=$(echo "$freqC" | wc -l );nFreqP0=$(echo "$freqP0" | wc -l );nFreqA=$(echo "$freqA" | wc -l );nFreqN=$(echo "$freqN" | wc -l )

if [[ $nFreqC -gt 1 && $nFreqP0 -gt 1  && $nFreqN -gt 1 && $nFreqA -gt 1 ]]; then
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tPopTE2\tND\tND\tND\tND" >> $DIR/ANALYSES/post_analysis/frequency_insertions_all_${class}_$type.tab

echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> $DIR/ANALYSES/post_analysis/fisher_PopTE2_all_${class}_$type.tab

else

if [[ $nFreqC -gt 1 ]];then
freqC=""
fi

if [[ $nFreqP0 -gt 1 ]];then
freqP0=""
fi

if [[ $nFreqN -gt 1 ]];then
freqN=""
fi

if [[ $nFreqA -gt 1 ]];then
freqA=""
fi


echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tPopTE2\t$freqP0\t$freqC\t$freqA\t$freqN" >> $DIR/ANALYSES/post_analysis/frequency_insertions_all_${class}_$type.tab

#### wilcox popte2
if [[ $type == "A" ]]; then
readsPopTE2A=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_A/${class}_A.teinsertions.family) | grep -E "($family.*){2}"| cut -f 13-14 | tr '/' '\t')
readsPopTE2N=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_N/${class}_N.teinsertions.family) | grep -E "($family.*){2}"| cut -f 13-14 | tr '/' '\t')
readsPopTE2C=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_C/${class}_C.teinsertions.family) | grep -E "($family.*){2}"| cut -f 13-14 | tr '/' '\t')
readsPopTE2P0=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_P0/${class}_P0.teinsertions.family) | grep -E "($family.*){2}"| cut -f 13-14 | tr '/' '\t'|head -n1)

nreadsPopTE2A=$(echo "$readsPopTE2A" | wc -l)
nreadsPopTE2N=$(echo "$readsPopTE2N" | wc -l)
nreadsPopTE2C=$(echo "$readsPopTE2C" | wc -l)
nreadsPopTE2P0=$(echo "$readsPopTE2P0" | wc -l)

if [[ $nreadsPopTE2A -gt 1 ]];then
readsPopTE2A=""
fi

if [[ $nreadsPopTE2N -gt 1 ]];then
readsPopTE2N=""
fi

if [[ $nreadsPopTE2C -gt 1 ]];then
readsPopTE2C=""
fi

if [[ $nreadsPopTE2P0 -gt 1 ]];then
readsPopTE2P0=""
fi

echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tA\tsel\t$readsPopTE2A" > $DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tN\tno_sel\t$readsPopTE2N" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tC\tC\t$readsPopTE2C" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tP0\tP0\t$readsPopTE2P0" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab

else
readsPopTE2A=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_A/${class}_A.teinsertions.family) | grep -E "($family.*){2}" | cut -f 13-14 | tr '/' '\t')
readsPopTE2N=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_N/${class}_N.teinsertions.family) | grep -E "($family.*){2}" | cut -f 13-14| tr '/' '\t')
readsPopTE2C=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_C/${class}_C.teinsertions.family) | grep -E "($family.*){2}" | cut -f 13-14| tr '/' '\t')
readsPopTE2P0=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(cut -f2- $DIR/ANALYSES/PoPoolationTE2/output/${class}_P0/${class}_P0.teinsertions.family) | grep -E "($family.*){2}" | cut -f 13-14| tr '/' '\t')

nreadsPopTE2A=$(echo "$readsPopTE2A" | wc -l)
nreadsPopTE2N=$(echo "$readsPopTE2N" | wc -l)
nreadsPopTE2C=$(echo "$readsPopTE2C" | wc -l)
nreadsPopTE2P0=$(echo "$readsPopTE2P0" | wc -l)

if [[ $nreadsPopTE2A -gt 1 ]];then
readsPopTE2A=""
fi

if [[ $nreadsPopTE2N -gt 1 ]];then
readsPopTE2N=""
fi

if [[ $nreadsPopTE2C -gt 1 ]];then
readsPopTE2C=""
fi

if [[ $nreadsPopTE2P0 -gt 1 ]];then
readsPopTE2P0=""
fi

echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tA\tno_sel\t$readsPopTE2A" > $DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tN\tsel\t$readsPopTE2N" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tC\tC\t$readsPopTE2C" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tP0\tP0\t$readsPopTE2P0" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab
fi

Rscript $DIR/ANALYSES/post_analysis/fisher_all.R "$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.PopTE2.tab" PopTE2 ${class} $type

fi
Rscript family_TEMP2.R $DIR/ANALYSES/TEMP2/output/${class}_C/insertion/${class}_C.insertion.bed
freqC=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_C/insertion/${class}_C.insertion.bed.family | cut -f 1-6,16) | grep -E "($family.*){2}" | cut -f 10)

Rscript family_TEMP2.R $DIR/ANALYSES/TEMP2/output/${class}_A/insertion/${class}_A.insertion.bed
freqA=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_A/insertion/${class}_A.insertion.bed.family | cut -f 1-6,16) | grep -E "($family.*){2}" | cut -f 10)

Rscript family_TEMP2.R $DIR/ANALYSES/TEMP2/output/${class}_N/insertion/${class}_N.insertion.bed
freqN=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_N/insertion/${class}_N.insertion.bed.family | cut -f 1-6,16) | grep -E "($family.*){2}" | cut -f 10)

Rscript family_TEMP2.R $DIR/ANALYSES/TEMP2/output/${class}_P0/insertion/${class}_P0.insertion.bed
freqP0=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_P0/insertion/${class}_P0.insertion.bed.family | cut -f 1-6,16) | grep -E "($family.*){2}" | cut -f 10)

nFreqC=$(echo "$freqC" | wc -l );nFreqP0=$(echo "$freqP0" | wc -l );nFreqA=$(echo "$freqA" | wc -l );nFreqN=$(echo "$freqN" | wc -l )

if [[ $nFreqC -gt 1 && $nFreqP0 -gt 1  && $nFreqN -gt 1 && $nFreqA -gt 1 ]]; then
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tTEMP2\tND\tND\tND\tND" >> $DIR/ANALYSES/post_analysis/frequency_insertions_all_${class}_$type.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> $DIR/ANALYSES/post_analysis/fisher_TEMP2_all_${class}_$type.tab
else

if [[ $nFreqC -gt 1 ]];then
freqC=""
fi

if [[ $nFreqP0 -gt 1 ]];then
freqP0=""
fi

if [[ $nFreqN -gt 1 ]];then
freqN=""
fi

if [[ $nFreqA -gt 1 ]];then
freqA=""
fi

echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tTEMP2\t$freqP0\t$freqC\t$freqA\t$freqN" >> $DIR/ANALYSES/post_analysis/frequency_insertions_all_${class}_$type.tab

### fisher temp2
if [[ $type == "A" ]]; then
readsTEMP2A=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_A/insertion/${class}_A.insertion.bed.family) | grep -E "($family.*){2}" | cut -f 10,13,14)
readsTEMP2N=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_N/insertion/${class}_N.insertion.bed.family) | grep -E "($family.*){2}" | cut -f 10,13,14)
readsTEMP2C=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_C/insertion/${class}_C.insertion.bed.family)| grep -E "($family.*){2}" | cut -f 10,13,14)
readsTEMP2P0=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_P0/insertion/${class}_P0.insertion.bed.family)| grep -E "($family.*){2}" | cut -f 10,13,14)

nreadsTEMP2A=$(echo "$readsTEMP2A" | wc -l)
nreadsTEMP2N=$(echo "$readsTEMP2N" | wc -l)
nreadsTEMP2C=$(echo "$readsTEMP2C" | wc -l)
nreadsTEMP2P0=$(echo "$readsTEMP2P0" | wc -l)

if [[ $nreadsTEMP2A -gt 1 ]];then
readsTEMP2A=""
fi

if [[ $nreadsTEMP2N -gt 1 ]];then
readsTEMP2N=""
fi

if [[ $nreadsTEMP2C -gt 1 ]];then
readsTEMP2C=""
fi

if [[ $nreadsTEMP2P0 -gt 1 ]];then
readsTEMP2P0=""
fi

echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tA\tsel\t$readsTEMP2A" > $DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tN\tno_sel\t$readsTEMP2N" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tC\tC\t$readsTEMP2C" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tP0\tP0\t$readsTEMP2P0" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab

else
readsTEMP2A=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_A/insertion/${class}_A.insertion.bed.family)| grep -E "($family.*){2}" | cut -f 10,13,14)
readsTEMP2N=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_N/insertion/${class}_N.insertion.bed.family)| grep -E "($family.*){2}" | cut -f 10,13,14)
readsTEMP2C=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_C/insertion/${class}_C.insertion.bed.family)| grep -E "($family.*){2}" | cut -f 10,13,14)
readsTEMP2P0=$(bedtools window -w 1 -a $DIR/ANALYSES/common_insertions/tmp/${class}_$type.$TE.$start.$end.all.bed -b <(tail -n +2 $DIR/ANALYSES/TEMP2/output/${class}_P0/insertion/${class}_P0.insertion.bed.family)| grep -E "($family.*){2}" | cut -f 10,13,14)

nreadsTEMP2A=$(echo "$readsTEMP2A" | wc -l)
nreadsTEMP2N=$(echo "$readsTEMP2N" | wc -l)
nreadsTEMP2C=$(echo "$readsTEMP2C" | wc -l)
nreadsTEMP2P0=$(echo "$readsTEMP2P0" | wc -l)

if [[ $nreadsTEMP2A -gt 1 ]];then
readsTEMP2A=""
fi

if [[ $nreadsTEMP2N -gt 1 ]];then
readsTEMP2N=""
fi

if [[ $nreadsTEMP2C -gt 1 ]];then
readsTEMP2C=""
fi

if [[ $nreadsTEMP2P0 -gt 1 ]];then
readsTEMP2P0=""
fi

echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tA\tno_sel\t$readsTEMP2A" > $DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tN\tsel\t$readsTEMP2N" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tC\tC\t$readsTEMP2C" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab
echo -e "${class}_$type\t$TE\t$chr\t$start\t$end\tP0\tP0\t$readsTEMP2P0" >>$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab
fi

Rscript $DIR/ANALYSES/post_analysis/fisher_all.R "$DIR/ANALYSES/post_analysis/tmp_reads/${class}_$type.$TE.$chr.$start.$end.reads.TEMP2.tab" TEMP2 ${class} $type

Rscript $DIR/ANALYSES/post_analysis/visualize_insertions_all.R "$TE" "$chr" "$start" "$end" "${class}_$type" "${class}_$type.$TE.$chr.$start.$end" ${class} $type

fi

#done <<< "$TEs"


# HEADER="Line\tCond\tTE\tchr\tstart\tend\tprogram\tP0\tC\tA\tN"
# sed -i "1 i$HEADER" $DIR/ANALYSES/post_analysis/frequency_insertions_all_${class}_$type.tab
# HEADER="Line\tTE\tchr\tstart\tend\tn\tp\tp.signif\tgroup1\tgroup2\tn\tp\tp.adj\tp.adj.signif"
# sed -i "1 i$HEADER" $DIR/ANALYSES/post_analysis/fisher_PopTE2_all_${class}_$type.tab
# sed -i "1 i$HEADER" $DIR/ANALYSES/post_analysis/fisher_TEMP2_all_${class}_$type.tab
# HEADER="Line\tTE\tchr\tstart\tend\tA\tN\tC\tP0\tFC"
# sed -i "1 i$HEADER" $DIR/ANALYSES/post_analysis/FC_PopTE2_all_${class}_$type.tab
# sed -i "1 i$HEADER" $DIR/ANALYSES/post_analysis/FC_TEMP2_all_${class}_$type.tab
