
rule kraken2: #
    input:
        se_fastq=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (wildcards.datatype,
                                                                                                    config["fastq_extension"])),
                    fileprefix=input_file_prefix_dict[wildcards.datatype],
                    allow_missing=True) if wildcards.datatype not in config["paired_fastq_based_data"] else [],
        forward_fastq=[], # TODO: implement contamination scan for hic and other types of pe data
        reverse_fastq=[], # TODO: implement contamination scan for hic and other types of pe data
        db=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["path"]
    output:
        summary=out_dir_path / "contamination_scan/kraken2/{datatype}/kraken2.{database}.report",
        out=out_dir_path / "contamination_scan/kraken2/{datatype}/kraken2.{database}.out",
    params:
        memory_mapping=lambda wildcards: "" if config["allowed_databases"]["kraken2"][wildcards.database]["in_memory"] else  " --memory-mapping ",
        paired=lambda wildcards: " --paired " if wildcards.datatype in config["paired_fastq_based_data"] else ""
    log:
        std=output_dict["log"]  / "kraken2.{database}.{datatype}.log",
        cluster_log=output_dict["cluster_log"] / "kraken2.{database}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "kraken2.{database}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "kraken2.{database}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["threads"] ,
        time=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["time"] ,
        mem=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["memory_mb"] ,
    threads: lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["threads"] ,

    shell:
        " kraken2 --threads {threads} {params.memory_mapping} {params.paired}  --db {input.db} "
        " --output {output.out} --report {output.summary} {input.se_fastq} > {log.std} 2>&1"

