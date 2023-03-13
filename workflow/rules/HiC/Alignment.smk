
rule bwa_index:
    input:
        fasta=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"])
    output:
        index=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta.bwt" % config["genome_name"])
    log:
        std=output_dict["log"]  / "bwa_index.purge_dups.{assembler}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "bwa_index.purge_dups.{assembler}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_index.purge_dups.{assembler}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_indexs.purge_dups.{assembler}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["bwa_index"] ,
        time=parameters["time"]["bwa_index"],
        mem=parameters["memory_mb"]["bwa_index"]
    threads: parameters["threads"]["bwa_index"]

    shell:
        " bwa index -a bwtsw {input.fasta} 1>{log.std} 2>&1;"

rule bwa_map: #
    input:
        index=rules.bwa_index.output.index,
        reference=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"]),
        fastq=output_dict["data"] / ("fastq/hic/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        bam=out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{fileprefix}.bam" % config["genome_name"])
    params:
        id="{0}_hic".format(config["genome_name"])
    log:
        map=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.sort.log",
        filter=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]

    shell:
        " bwa mem -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:Illumina\\tLB:x\' "
        " {input.reference} {input.fastq} 2>{log.map} | filter_five_end.pl 2>{log.filter} | samtools view -Sb - > {output.bam} 2>{log.sort} "

"""
rule bam_merge: # TODO: add nanopore support
    input:
        index=rules.bwa_index.output.index,
        reference=out_dir_path  / ("{assembly_stage}/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"]),
        fastq=output_dict["data"] / ("fastq/hic/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        bam=out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{fileprefix}.bam" % config["genome_name"])
    params:
        id="{0}_hic".format(config["genome_name"])
    log:
        map=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]

    shell:
        " bwa mem -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:Illumina\\tLB:x\' "
        " {input.reference} {input.fastq} 2>{log.map} | perl filter_five_end.pl | samtools view -Sb - > {output.bam} 2>{log.sort} "



echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/$SRA\_1.bam $FILT_DIR/$SRA\_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA.bam -



#! /bin/bash

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

SRA='basename_of_fastq_files'
LABEL='overall_exp_name'
BWA='/software/bwa/bwa-0.7.12/bwa'
SAMTOOLS='/software/samtools/samtools-1.3.1/samtools'
IN_DIR='/path/to/gzipped/fastq/files'
REF='/path/to/reference_sequences/reference_sequeneces.fa'
FAIDX='$REF.fai'
PREFIX='bwa_index_name'
RAW_DIR='/path/to/write/out/bams'
FILT_DIR='/path/to/write/out/filtered/bams'
FILTER='/path/to/filter_five_end.pl'
COMBINER='/path/to/two_read_bam_combiner.pl'
STATS='/path/to/get_stats.pl'
PICARD='/software/picard/picard-2.6.0/build/libs/picard.jar'
TMP_DIR='/path/to/write/out/temporary/files'
PAIR_DIR='/path/to/write/out/paired/bams'
REP_DIR='/path/to/where/you/want/deduplicated/files'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='/path/to/final/merged/alignments/from/any/biological/replicates'
MAPQ_FILTER=10
CPU=12

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
$BWA index -a bwtsw -p $PREFIX $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
$BWA mem -t $CPU $REF $IN_DIR/$SRA\_1.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
$BWA mem -t $CPU $REF $IN_DIR/$SRA\_2.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_2.bam


echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

"""
