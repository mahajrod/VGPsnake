
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

rule ref_faidx:
    input:
        reference=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"]),
    output:
        fai=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fai" % config["genome_name"]),
    log:
        std=output_dict["log"]  / "ref_faidx.purge_dups.{assembler}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "ref_faidx.purge_dups.{assembler}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ref_faidx.purge_dups.{assembler}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ref_faidx.purge_dups.{assembler}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["ref_faidx"] ,
        time=parameters["time"]["ref_faidx"],
        mem=parameters["memory_mb"]["ref_faidx"]
    threads: parameters["threads"]["ref_faidx"]

    shell:
        " samtools faidx -o {output.fai} {input.reference} 1>{log.std} 2>&1;"

rule ref_dict:
    input:
        reference=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"]),
    output:
        dict=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.dict" % config["genome_name"]),
    log:
        std=output_dict["log"]  / "ref_dict.purge_dups.{assembler}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "ref_dict.purge_dups.{assembler}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ref_dict.purge_dups.{assembler}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ref_dict.{assembler}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["ref_dict"] ,
        time=parameters["time"]["ref_dict"],
        mem=parameters["memory_mb"]["ref_dict"]
    threads: parameters["threads"]["ref_dict"]

    shell:
        " picard CreateSequenceDictionary R={input.reference} O={output.dict} > {log.std} 2>&1;"

"""
rule index_bam:
    input:
        bam="{bam_path,.*}"
    output:
        "{bam_path,.*}.bai"
    log:
        std=output_dict["log"]  / "ref_dict.purge_dups.{assembler}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "ref_dict.purge_dups.{assembler}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ref_dict.purge_dups.{assembler}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ref_dict.{assembler}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["ref_dict"] ,
        time=parameters["time"]["ref_dict"],
        mem=parameters["memory_mb"]["ref_dict"]
    threads: parameters["threads"]["ref_dict"]
    shell:
        "samtools index -@ {threads} {input} > {log.std} 2>&1"
"""