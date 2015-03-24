#!/bin/bash

### Shell scripted version of exome analysis pipeline for trios.
### GATK Genotype refinement protocol was referred to during the 
### development of this workflow.
### http://gatkforums.broadinstitute.org/discussion/4723/genotype-refinement-workflow
### This script will be submitted three times to Condor.  Subsequently, we will gather
### the output of each submission and call variants based on all three samples with pedigree info.
### bash trio_pipeline.sh <mate1> <mate2> <sample_name> <child, mother, or father>
### John Letaw 03/10/15

### Specify number of worker threads.
DTHREADS=20  ### Specify this.
CTHREADS=8

### Command line arguments
MATE1=$1
MATE2=$2
SAM_BASE=$(basename $MATE1 .fastq)
OUT_SAM=$SAM_BASE.sam
SAMPLE=$3
REL=$4

### Program directories
BIN_DIR="/opt/installed"
PICARD="$BIN_DIR/picard/picard-tools-1.110"
GATK="$BIN_DIR/GATK/GenomeAnalysisTK-3.1.jar"

### Project directories
PROJ="/home/groups/clinical/RichardsLab/letaw/trio_pipeline"  ### Specify this.
RESOURCE="/home/groups/clinical/RichardsLab/resources"
RAW="$PROJ/raw"
LOGS="$PROJ/logs"
BWA="$PROJ/bwa"
SORTED="$PROJ/sorted"
DUPS="$PROJ/mark_dups"
INDEL="$PROJ/indel_realigned"
BQSR="$PROJ/bqsr"
HC="$PROJ/hc"
GG="$PROJ/genotype_gvcfs"
VR="$PROJ/var_recal"
PHASED="$PROJ/phased"
DOC="$PROJ/doc"

### Reference files
HG="$RESOURCE/genomes/hg19/ucsc.hg19.fasta"
#MILLS="/home/users/letaw/clinical/RichardsLab/resources/variation/mills/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
MILLS="$RESOURCE/variation/mills/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
#DBSNP="/home/users/letaw/clinical/RichardsLab/resources/variation/dbsnp/142/grch37_p13_142.vcf"
DBSNP="$RESOURCE/variation/dbsnp/138/dbsnp_138.hg19.vcf"
THG="$RESOURCE/variation/1000g/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
HAPMAP="$RESOURCE/variation/hapmap/hapmap_3.3.hg19.sites.vcf"
OMNI="$RESOURCE/variation/omni/1000G_omni2.5.hg19.sites.vcf"
INTERVALS="/home/users/letaw/clinical/RichardsLab/resources/intervals/WES.interval_list"

### NA12878
#@NS500390:8:H2CWJAFXX:1:11101:2206:985 1:N:0:ACAGTG
### NA12891
#@NS500390:8:H2CWJAFXX:1:11101:17980:985 1:N:0:GCCAAT
### NA12892
#@NS500390:8:H2CWJAFXX:1:11101:3349:985 1:N:0:CAGATC

### Read group information.  This will be auto-created from FASTQ ID lines.
RGID=$(head -1 $MATE1 | cut -d ':' -f 1-4)
RGSM=$SAMPLE
RGPL="illumina"  ### Hard code since we will be only dealing with Illumina data through this pipeline.
RGLB=$(head -1 $MATE1 | cut -d ':' -f 3-4)
RGPU=$(head -1 $MATE1 | cut -d ':' -f 10)
RGDS=$REL
RG="@RG\tID:$RGID\tSM:$RGSM\tPL:$RGPL\tLB:$RGLB\tPU:$RGPU\tDS:$RGDS"

### Create directory structure.
mkdir -p $RAW $BWA $SORTED $DUPS $INDEL $BQSR $HC $GG $VR $PHASED $DOC

### Index your genome if you haven't done so already with bwa index.
if [ ! -f $HG.sa ]
then
    echo "Indexing $HG."
    bwa index $HG
else
    echo "$HG is already indexed."
fi

# ### Align to reference.
# echo "Aligning to $HG."
# $BIN_DIR/bwa mem -M -R $RG -t $THREADS $HG $MATE1 $MATE2 > $BWA/$OUT_SAM

# ### Sort your SAM file.  Do not produce an index at this stage.
# echo "Sorting the SAM and outputting a BAM."
# java -jar -Xmx64g $PICARD/SortSam.jar MAX_RECORDS_IN_RAM=10000000 SO=coordinate I=$BWA/$OUT_SAM O=$SORTED/$SAM_BASE\_sorted.bam

# ### Mark duplicates, create and index.
# echo "Marking duplicates."
# java -jar -Xmx64g $PICARD/MarkDuplicates.jar MAX_RECORDS_IN_RAM=10000000 CREATE_INDEX=TRUE I=$SORTED/$SAM_BASE\_sorted.bam O=$DUPS/$SAM_BASE\_sorted_dups.bam M=$DUPS/$SAM_BASE.metrics

# ### Create the interval file, if it doesn't already exist.
# if [ ! -f $PROJ/target_intervals.list ] && [ "$4" == "child" ] 
# then
#     echo "Creating intervals for IndelRealigner."
#     java -jar $GATK -R $HG -T RealignerTargetCreator -I $DUPS/$SAM_BASE\_sorted_dups.bam -known $MILLS -o $PROJ/target_intervals.list -nt $DTHREADS
# else
#     while [ ! -f $PROJ/target_intervals.list ]
#     do
# 	sleep 5
#     done
#     echo "Interval file has already been created, proceeding to IndelRealigner."
# fi

# ### Run a realignment around targeted indels.
# echo "Realigning around targeted indel regions."
# java -jar $GATK -R $HG -T IndelRealigner -I $DUPS/$SAM_BASE\_sorted_dups.bam -targetIntervals $PROJ/target_intervals.list -known $MILLS -o $INDEL/$SAM_BASE\_sorted_dups_realigned.bam

# ### Base recalibration.
# echo "Base recalibration round 1."
# java -jar $GATK -R $HG -T BaseRecalibrator -I $INDEL/$SAM_BASE\_sorted_dups_realigned.bam -knownSites $DBSNP -knownSites $MILLS -knownSites $THG -o $BQSR/$SAM_BASE\_recal_data.table -nct $CTHREADS

# echo "Base recalibration round 2."
# java -jar $GATK -R $HG -T BaseRecalibrator -I $INDEL/$SAM_BASE\_sorted_dups_realigned.bam -knownSites $DBSNP -knownSites $MILLS -knownSites $THG -BQSR $BQSR/$SAM_BASE\_recal_data.table -o $BQSR/$SAM_BASE\_post_recal_data.table -nct $CTHREADS

# echo "Base recalibration plots."
# java -jar $GATK -R $HG -T AnalyzeCovariates -before $BQSR/$SAM_BASE\_recal_data.table -after $BQSR/$SAM_BASE\_post_recal_data.table -plots $BQSR/$SAM_BASE\_recal_plot.pdf

# echo "Applying recalibration."
# java -jar $GATK -R $HG -T PrintReads -I $INDEL/$SAM_BASE\_sorted_dups_realigned.bam -BQSR $BQSR/$SAM_BASE\_post_recal_data.table -o $BQSR/$SAM_BASE\_sorted_dups_realigned_bqsr.bam -nct $CTHREADS

# ### Calling variants.
# echo "Calling variants with HaplotypeCaller."
# java -jar $GATK -R $HG -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -I $BQSR/$SAM_BASE\_sorted_dups_realigned_bqsr.bam -o $HC/$SAM_BASE.gvcf -nct $CTHREADS

# ### Merge and genotype gVCFs.

# if [ "$4" == "child" ]
# then
#     echo "Merging files."

#     for log in $(ls -1 $LOGS/*.{1,2}.log)
#     do
# 	echo "Waiting for $log."
# 	until grep -q "return value 100" $log
# 	do
# 	    sleep 5
# 	done
#     done

#     java -jar $GATK -R $HG -T GenotypeGVCFs `for line in $(ls -1 $HC/*.gvcf); do echo -n "-V $line "; done` -D $DBSNP -o $GG/trio.vcf -nt $DTHREADS
# else
#     echo "DONE"
#     exit 100
# fi

# ### Build SNP recalibration model.
# echo "Building the SNP recalibration model for VQSR."
# java -jar $GATK -R $HG -T VariantRecalibrator -input $GG/trio.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $THG -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an DP -an QD -an FS -an MQ -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $VR/recalibrate_SNP.recal -tranchesFile $VR/recalibrate_SNP.tranches -rscriptFile $VR/recalibrate_SNP_plots.R -nt $DTHREADS

# ### Apply recalibration.
# echo "Applying VQSR."
# java -jar $GATK -R $HG -T ApplyRecalibration -input $GG/trio.vcf -mode SNP --ts_filter_level 99.0 -recalFile $VR/recalibrate_SNP.recal -tranchesFile $VR/recalibrate_SNP.tranches -o $VR/trio_recal_snps_raw_indels.vcf -nt $DTHREADS

# ### Build INDEL recalibration model
# echo "Building the indel recalibration model for VQSR."
# java -jar $GATK -T VariantRecalibrator -R $HG -input $VR/trio_recal_snps_raw_indels.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS -an QD -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $VR/recalibrate_INDEL.recal -tranchesFile $VR/recalibrate_INDEL.tranches -rscriptFile $VR/recalibrate_INDEL_plots.R -nt $DTHREADS 

# ### Apply recalibration.
# echo "Applying VQSR."
# java -jar $GATK -T ApplyRecalibration -R $HG -input $VR/trio_recal_snps_raw_indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile $VR/recalibrate_INDEL.recal -tranchesFile $VR/recalibrate_INDEL.tranches -o $VR/trio_recalibrated_variants.vcf 

# ### Apply phasing based on pedigree info.
# echo "Applying phasing."
# java -jar $GATK -T PhaseByTransmission -R $HG -V $VR/trio_recalibrated_variants.vcf -ped $PROJ/raw/pedigree -o $PHASED/final_trio.vcf -mvf $VR/trio_mend_viol

### Depth Of Coverage.
echo "Running DepthOfCoverage."
java -jar $GATK -R $HG -T DepthOfCoverage -I $BQSR/$SAM_BASE\_sorted_dups_realigned_bqsr.bam -pt sample --outputFormat rtable -L $INTERVALS -isr UNION -baqGOP 40 -DBQ 0 -S STRICT -im ALL --maxBaseQuality 127 -mbq 20 --maxMappingQuality 2147483647 -mmq 20 --nBins 499 -omitIntervals --start 1 --stop 500 -o $DOC/$SAM_BASE\_doc
