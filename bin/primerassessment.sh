#!/bin/bash

version="0.20230630"
USAGE="
primerassessment.sh is meant to extract paired-end reads in the coordinates specified by a bed file of amplicon start and end coordinates.

Usage:
primerassessment.sh
"
############################################################

bam=test/SRR13957125.sorted.bam
bed=/home/eriny/sandbox/Cecret/schema/artic_V3_nCoV-2019.insert.bed 

i=0


header="sample"
result="test"

while read line
do
    ref=$(echo   $line | awk '{print $1}' )
    start=$(echo $line | awk '{print $2}' )
    stop=$(echo  $line | awk '{print $3}' )
    name=$(echo  $line | awk '{print $4}' )
    remove_bed=remove_${i}.bed
    header="$header,$name"

    echo "working on line $i"
    echo "determining coverage for $line with ref = $ref, start = $start, stop = $stop, and name = $name"

    echo -e "$ref\t1\t$((start - 1))"      >  $remove_bed
    echo -e "$ref\t$((stop + 1))\t1000000" >> $remove_bed

    # samtools view -b test/SRR13957125.sorted.bam -f2 MN908947.3:47-447 -h | samtools view - -U MN908947.3:1-46 -h | samtools view - -U MN908947.3:448-400000 -h | head

    # remove fastq files that map outside of the region
    echo "removing reads that map outside of region"
    samtools view -bh $bam -U test/test_${i}_remove.bam -o test/test_${i}_inregion1.bam -L $remove_bed
    samtools index test/test_${i}_remove.bam

    echo "viewing reads in region"
    samtools view -bh -f2 test/test_${i}_remove.bam $ref:$start-$stop -o test/test_${i}_mapped.bam
    samtools index test/test_${i}_mapped.bam

    echo "getting meandepth"
    meandepth=$(samtools coverage test/test_${i}_mapped.bam -r $ref:$start-$stop | tail -n 1 | awk '{print $7}' )

    result="$result,$meandepth"
    echo $header
    echo $result

    i=$((i + 1))

done < $bed

echo $header >  primer_assessment.csv
echo $result >> primer_assessment.csv