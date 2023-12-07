#!/bin/bash

# Run by using sh align_single_fastq_to_reference.sh FASTQ_FILE_R1_001.fastq.gz REFERENCE_FILE.fasta

FILE=$1
REF=$2

PREFIX=${FILE%_R1_001.fastq.gz}

echo "File prefix: " $PREFIX
echo "Reference file: " $REF

# Copy reference file over... Kinda janky.
# cp $REF .

##################### The code starts here
# Base name of the reference file.
REF_BASENAME=$(basename $REF)
REF_BASENAME=${REF_BASENAME%.fa*}

REF_INDEX_NAME=${REF%.fa*}

FOR_READ=${PREFIX}_R1_001.fastq.gz
REV_READ=${PREFIX}_R2_001.fastq.gz

# Build an index for the reference. This only needs to be done once.
echo -e "########################\n${PREFIX}: BUILDING BOWTIE INDEX FOR REFERENCE \n########################"
bowtie2-build ${REF} ${REF_INDEX_NAME}

# First trim the reads using BBDuk.
echo -e "########################\n${PREFIX}: TRIM USING BBDUK \n########################"
bbduk.sh in=${PREFIX}_R1_001.fastq.gz in2=${PREFIX}_R2_001.fastq.gz \
    out=${PREFIX}_R1_clean_TMP.fastq.gz out2=${PREFIX}_R2_clean_TMP.fastq.gz \
    qtrim=rl trimq=28 \
    maq=25 \
    overwrite=true

# Convert bases below q28 to be N.
echo -e "########################\n${PREFIX}: MASK TO Ns \n########################"
/broad/thechenlab/Dawn/software/seqtk/seqtk seq -q28 -n N ${PREFIX}_R1_clean_TMP.fastq.gz > ${PREFIX}_R1_clean.fastq.gz
/broad/thechenlab/Dawn/software/seqtk/seqtk seq -q28 -n N ${PREFIX}_R2_clean_TMP.fastq.gz > ${PREFIX}_R2_clean.fastq.gz

echo -e "########################\n${PREFIX}: ALIGN TO REFERENCE \n########################"
bowtie2 -t -p 8 -x ${REF_INDEX_NAME} \
    -1 ${PREFIX}_R1_clean.fastq.gz \
    -2 ${PREFIX}_R2_clean.fastq.gz \
    --local \
    --minins 0 --maxins 2500 \
    -S ${PREFIX}_aligned_to_${REF_BASENAME}.sam

echo -e "########################\n${PREFIX}: CONVERT TO BAM \n########################"
# Need to sort and build index. Required for downstream.
# samtools view -bS ${PREFIX}_aligned_to_${REF_BASENAME}.sam >${PREFIX}_aligned_to_${REF_BASENAME}.bam
samtools sort -@ 4 ${PREFIX}_aligned_to_${REF_BASENAME}.sam -o ${PREFIX}_aligned_to_${REF_BASENAME}.bam -O bam

# Remove SAM file because it's big.
rm ${PREFIX}_aligned_to_${REF_BASENAME}.sam
# Also remove the temp fastqs.
rm ${PREFIX}_R1_clean_TMP.fastq.gz
rm ${PREFIX}_R2_clean_TMP.fastq.gz
rm ${PREFIX}_R1_clean.fastq.gz
rm ${PREFIX}_R2_clean.fastq.gz

# Index the bam file.
samtools index ${PREFIX}_aligned_to_${REF_BASENAME}.bam

# Make pileup table.
echo -e "########################\n${PREFIX}: MAKE PILEUP TABLE \n########################"
python get_pileup_table_hace_amplicon.py -b ${PREFIX}_aligned_to_${REF_BASENAME}.bam -r $REF

echo -e "########################\n${PREFIX}: DONE! \n########################"

