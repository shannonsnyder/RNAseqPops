#!/bin/bash                                                                                       

#SBATCH -n 1                        # number of cores                                            
#SBATCH -t 0-12:00                  # wall time (D-HH:MM)                                        
##SBATCH -A ssnyde11                # Account hours will be pulled from (commented out with double # in front)                                                                                    
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)                                        
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)                                        
#SBATCH --mail-type=ALL             # Send a notification



BWAdir=/home/ssnyde11/scratch/Daphnia_RNAseq_090719
GENOME_DIR=/home/ssnyde11/DpGENOME
GENOME_FILE=PA42.4.1.fasta
TAGDUST=/home/ssnyde11/genome_analysis/tagdust-2.33/src/tagdust
BWA=/packages/7x/bwakit/0.7.12/bwa.kit/bwa
RNAfile=/home/ssnyde11/DpGENOME/Dpul-rDNA.fa
fastqDir=/home/ssnyde11/scratch/Daphnia_RNAseq_090719
SAMTOOLS=/packages/7x/samtools/1.8/bin/samtools

Removing ribosomal RNA hits with tagdust:


Removing ribosomal RNA hits with tagdust:



THREADS=8

SAMPLE1=GSF2349-LPB-2014-32-5_S5
SAMPLE2=GSF2349-LPB-2014-32-2-new_S6
SAMPLE3=GSF2349-LPA-2014-32-5-new_S4
SAMPLE4=GSF2349-LPA-2014-16-REP2_S2
SAMPLE5=GSF2349-LPA-2014-16-REP1_S1
SAMPLE6=GSF2349-LPA-2014-16-3_S3
SAMPLE7=GSF2349-POV-2014-12-REP4_S12
SAMPLE8=GSF2349-POV-2014-12-REP3_S11
SAMPLE9=GSF2349-POV-2014-12-REP2_S10
SAMPLE10=GSF2349-POV-2014-12-REP1_S9
SAMPLE11=GSF2349-NFL3-REP4_S20
SAMPLE12=GSF2349-NFL3-REP3_S19
SAMPLE13=GSF2349-NFL3-REP2_S18
SAMPLE14=GSF2349-NFL3-REP1_S17
SAMPLE15=GSF2349-LPB-2014-32-4-old_S8
SAMPLE16=GSF2349-LPB-2014-32-3_S7
SAMPLE17=GSF2349-KAP-2014-114-4-new_S14
SAMPLE18=GSF2349-KAP-2013-114-REP1_S13
SAMPLE19=GSF2349-KAP-2013-114-4_S16
SAMPLE20=GSF2349-KAP-2013-114-3_S15

daphnia_rep1_R1=${SAMPLE1}_R1_001.fastq
#daphnia_rep1_R2_${SAMPLE1}_R2_001.fastq
daphnia_rep2_R1=${SAMPLE2}_R1_001.fastq
#daphnia_rep2_R2=GSF2251-RTO_${SAMPLE2}_R2_001.fastq
daphnia_rep3_R1=${SAMPLE3}_R1_001.fastq
#daphnia_rep3_R2=GSF2251-RTO_${SAMPLE3}_R2_001.fastq
daphnia_rep4_R1=${SAMPLE4}_R1_001.fastq
#daphnia_rep4_R2=GSF2251-RTO_${SAMPLE4}_R2_001.fastq
daphnia_rep5_R1=${SAMPLE5}_R1_001.fastq
#daphnia_rep5_R2=GSF2251-RTO_${SAMPLE5}_R2_001.fastq
daphnia_rep6_R1=${SAMPLE6}_R1_001.fastq
#daphnia_rep6_R2=GSF2251-RTO_${SAMPLE6}_R2_001.fastq
daphnia_rep7_R1=${SAMPLE7}_R1_001.fastq
#daphnia_rep7_R2=GSF2251-RTO_${SAMPLE7}_R2_001.fastq
daphnia_rep8_R1=${SAMPLE8}_R1_001.fastq
#daphnia_rep8_R2=GSF2251-RTO_${SAMPLE8}_R2_001.fastq
daphnia_rep9_R1=${SAMPLE9}_R1_001.fastq
daphnia_rep10_R1=${SAMPLE10}_R1_001.fastq
daphnia_rep11_R1=${SAMPLE11}_R1_001.fastq
daphnia_rep12_R1=${SAMPLE12}_R1_001.fastq
daphnia_rep13_R1=${SAMPLE13}_R1_001.fastq
daphnia_rep14_R1=${SAMPLE14}_R1_001.fastq
daphnia_rep15_R1=${SAMPLE15}_R1_001.fastq
daphnia_rep16_R1=${SAMPLE16}_R1_001.fastq
daphnia_rep17_R1=${SAMPLE17}_R1_001.fastq
daphnia_rep18_R1=${SAMPLE18}_R1_001.fastq
daphnia_rep19_R1=${SAMPLE19}_R1_001.fastq
daphnia_rep20_R1=${SAMPLE20}_R1_001.fastq


OP1=${SAMPLE1}_trno_tagdusted
OP2=${SAMPLE2}_trno_tagdusted
OP3=${SAMPLE3}_trno_tagdusted
OP4=${SAMPLE4}_trno_tagdusted

OP5=${SAMPLE5}_trno_tagdusted
OP6=${SAMPLE6}_trno_tagdusted
OP7=${SAMPLE7}_trno_tagdusted
OP8=${SAMPLE8}_trno_tagdusted

OP9=${SAMPLE9}_trno_tagdusted
OP10=${SAMPLE10}_trno_tagdusted
OP11=${SAMPLE11}_trno_tagdusted
OP12=${SAMPLE12}_trno_tagdusted

OP13=${SAMPLE13}_trno_tagdusted
OP14=${SAMPLE14}_trno_tagdusted
OP15=${SAMPLE15}_trno_tagdusted
OP16=${SAMPLE16}_trno_tagdusted

OP17=${SAMPLE17}_trno_tagdusted
OP18=${SAMPLE18}_trno_tagdusted
OP19=${SAMPLE19}_trno_tagdusted
OP20=${SAMPLE20}_trno_tagdusted

cd $fastqDir

${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep1_R1}  -o ${OP1}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep2_R1} -o ${OP2}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep3_R1}  -o ${OP3}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep4_R1}  -o ${OP4}

${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep5_R1}  -o ${OP5}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep6_R1}  -o ${OP6}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep7_R1}  -o ${OP7}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep8_R1}  -o ${OP8}

${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep9_R1}  -o ${OP9}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep10_R1}  -o ${OP10}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep11_R1}  -o ${OP11}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep12_R1}  -o ${OP12}

${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep13_R1}  -o ${OP13}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep14_R1}  -o ${OP14}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep15_R1}  -o ${OP15}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep16_R1}  -o ${OP16}

${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep17_R1}  -o ${OP17}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep18_R1}  -o ${OP18}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep19_R1}  -o ${OP19}
${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${daphnia_rep20_R1}  -o ${OP20}


cd /home/ssnyde11/scratch/Daphnia_RNAseq_090719

echo "Indexing the Daphnia genome using bwa ..."
echo "bwa index ${GENOME_DIR}/${GENOME_FILE}"
${BWA}  index ${GENOME_DIR}/${GENOME_FILE}

echo "Starting alignments ..."
for fq in *_trno_tagdusted.fq;
do

        echo "bwa aln -t ${THREADS} -n 3 ${GENOME_DIR}/${GENOME_FILE} -f $(basename ${fq} _trno_tagdusted_READ1.fq).sai ${fq} $(basename ${fq} _READ1.fq)_READ2.fq"
  	${BWA} aln -t ${THREADS} -n 3 ${GENOME_DIR}/${GENOME_FILE} -f $(basename ${fq} _trno_tagdusted.fq).sai ${fq} 

	echo "$fq"

echo "${BWA} samse ${GENOME_DIR}/${GENOME_FILE} $(basename $fq .fq_trno_tagdusted.fq).sai $fq | ${SAMTOOLS} view -uS - | ${SAMTOOLS} sort -O BAM - > $(basename $fq .fq_trno_tagdusted.fq)_sorted.bam"
${BWA} samse ${GENOME_DIR}/${GENOME_FILE} $(basename $fq _trno_tagdusted.fq).sai $fq | ${SAMTOOLS} view -uS - | ${SAMTOOLS} sort -O BAM - > $(basename $fq _trno_tagdusted.fq)_sorted.bam

echo "samtools index -b $(basename $fq .fq_trno_tagdusted.fq)_sorted.bam "
${SAMTOOLS} index -b $(basename $fq .fq_trno_tagdusted.fq)_sorted.bam

#FILTERED_BAM=$(basename $fq _trno_tagdusted.fq)_filtered.bam
#.. post-alignment filtering for proper alignments and MAPQ >= 10:
#
echo "${SAMTOOLS} view -f 2 -q 10 -u ${SORTED_BAM} | ${SAMTOOLS} sort -O BAM -@ 10 - > $(basename $fq _trno_tagdusted.fq)_filtered.bam"
${SAMTOOLS} view -h -F 4 -q 10 -u $(basename $fq _trno_tagdusted.fq)_sorted.bam | ${SAMTOOLS} sort -O BAM -@ 8 - > $(basename $fq _trno_tagdusted.fq)_filtered.bam

echo "samtools index -b $(basename $fq .fq_trno_tagdusted.fq)_filtered.bam"
${SAMTOOLS} index -b $(basename $fq _trno_tagdusted.fq)_filtered.bam

done
