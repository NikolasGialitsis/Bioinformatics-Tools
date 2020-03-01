#!/bin/bash
#@AUTHOR GIALITSIS NIKOLAOS
echo "========= PROGRAM START ============="
echo "=============================       1          =============================="
r1=$(grep -o '@' sample_1.fastq | wc -l) #count number of reads in sample_1
r2=$(grep -o '@' sample_2.fastq | wc -l) #count number of reads in sample_2
echo "Number of reads in sample_1 is ${r1} and number of reads in sample_2 is ${r2}"

csam=$(samtools view -c -S aligned.sam)
echo "=============================       2          =============================="
echo "Number of reads in aligned.sam file is ${csam}"
mapped=$(samtools view -S -f 0x4 aligned.sam | wc -l)
echo "Number of mapped reads in aligned.sam file is ${mapped}"
unmapped=$(samtools view -S -F 0x4 aligned.sam | wc -l)
echo "Number of unmapped reads in aligned.sam file is ${unmapped}"
echo "=============================       3          =============================="


BAM=sample.bam
echo "Convert aligned.sam to ${BAM} file..."
samtools view -S -b aligned.sam > ${BAM}
if [ -f $BAM ]; then
    echo -e "\t succesfully created ${BAM}" 
else 
    echo -e "\t failed to create {$BAM}"
fi
echo "=============================       4          =============================="

SORTED=sorted
echo "Sorting ${BAM}..."
samtools sort ${BAM} ${SORTED} 
SORTED=sorted.bam
if [ -f ${SORTED} ]; then
    echo -e "\t succesfully created ${SORTED}" 
else 
    echo -e "\t failed to create {$SORTED}"
fi

size1=$(stat -c%s "${BAM}")
size2=$(stat -c%s "${SORTED}")

csort=$(samtools view -c ${SORTED})
echo "Number of reads in ${SORTED} file is ${csort}"

if [ ${size1} != ${size2} ];then
	echo -e "\t sizes of ${BAM} and ${SORTED} are different"
	if [ ${size1} -gt ${size2} ]; then
		echo -e "\t\t> ${BAM} is larger than ${SORTED} because a sorted bam file\n\t\t> can be further compressed than an unsorted one"
		if [ ${csam} -eq ${csort} ];then
			echo -e "\t\t\t>> however the number of reads in the files is the same\n\t\t\t>> so no information is lost"
		fi
	else
		echo -e "\t\t>> ${SORTED} is larger than ${BAM}"
	fi
fi

echo "=============================       5          =============================="
INDEX=${SORTED}.bai
echo "Indexing ${SORTED} to produce ${INDEX}"
samtools index ${SORTED}
if [ -f ${INDEX} ]; then
    echo -e "\t succesfully created ${INDEX}" 
else 
    echo -e "\t failed to create {$INDEX}"
fi
BAMIND=${BAM}.bai
echo "Indexing ${BAM} to produce ${BAMIND}"
samtools index ${BAM} 2> /dev/null 
if [ -f ${BAMIND} ]; then
    echo -e "\t succesfully created ${BAMIND}" 
else  
    echo -e "\t == failed to create ${BAMIND} ==\n\t because it is only possible to index BAM files on\n\t position, and only when the data is initially sorted by position "
fi

echo "=============================       6          =============================="
GENOME=human_g1k_v37.genome
echo -e "Calculating coverage of ${SORTED} on ${GENOME}..."
bedtools genomecov -ibam ${SORTED} -g ${GENOME} > genome_coverage
echo -e "...done\n\t >> output stored in genome_coverage file"

echo "=============================       7          =============================="
echo -e "Calculating coverage of ${SORTED} on ${GENOME} [BEDGRAPH FORMAT]..."
bedtools genomecov -bg -ibam ${SORTED} -g ${GENOME} > genome_coverage_bedgraph
echo -e "...done\n\t >> output stored in genome_coverage_bedgraph file"


echo "=============================       8          =============================="
FORWARD=forward_strand.out
REVERSE=reverse_strand.out
echo "Split ${BAM} into forward and reverse strands..."
echo -e "\t retrieve forward strand"
samtools view -F 0x10 ${BAM} -o ${FORWARD}
echo -e "\t\tforward strand stored in file \"${FORWARD}\""

echo -e "\t retrieve reverse strand"
samtools view -f 0x10 ${BAM} -o ${REVERSE}
echo -e "\t\treverse strand stored in file \"${REVERSE}\""


echo "=============================       9         =============================="

TARGET=TargetRegion.bed
TARGET100=TargetRegion.100bp.bed 
echo "Add 100bp at the end of ${TARGET} and save results to ${TARGET100}"
bedtools slop -i  ${TARGET} -g ${GENOME}  -l 0 -r 100 > ${TARGET100}
echo "=============================       10         =============================="

MERGED=TargetRegion.100bp.merged.bed
echo -e "Sorting by chromosome and then feature size in ascending order\n and merging ${TARGET100} "
bedtools sort -i ${TARGET100} > temp
mv temp ${TARGET100}
echo -e "\t> sorted"
bedtools merge -i ${TARGET100} > ${MERGED}
echo -e "\t> merged"


echo "=============================       11        =============================="
echo "Intersect ${SORTED} and ${MERGED}..."
INTERSECT=intersect.out
bedtools intersect -abam ${SORTED} -b ${MERGED} > ${INTERSECT}
echo "...done"
echo "There exist $(cat ${INTERSECT} | wc -l)  reads within ${INTERSECT}"
echo "=============================       12        =============================="
	

A=A.vcf 
B=B.vcf 
echo "Compress ${A} and ${B}..."
bgzip ${A} 2> /dev/null
bgzip ${B} 2> /dev/null
echo "...done"
echo "Index ${A} and ${B}..."
A=${A}.gz
B=${B}.gz
tabix -f -p vcf ${A} 
tabix -f -p vcf ${B}
echo "...done"

echo "Compare files ${A} and ${B} using vcf-compare..."
vcf-compare ${A} ${B} > compare.out
echo "...done"
echo -e "\t comparison results stored in \"compare.out\""

echo "========= PROGRAM FINISH ============="
exit 0