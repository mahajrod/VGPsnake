localrules: create_fastq_links

rule create_fastq_links:
    priority: 1000
    input:
        directory(input_dir_path / "{datatype}/fastq")

    output:
        directory(output_dict["data"] / "{datatype}/fastq/")
    log:
        std=output_dict["log"] / "preprocessing/{datatype}/create_fastq_links.log",
        cluster_log=output_dict["cluster_log"] / "preprocessing/{datatype}/create_fastq_links.cluster.log",
        cluster_err=output_dict["cluster_error"] / "preprocessing/{datatype}/create_fastq_links.cluster.err",
    benchmark:
        output_dict["benchmark"] / "preprocessing/{datatype}/create_fastq_links.benchmark.txt",
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["create_fastq_links"],
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"],
    threads:
        parameters["threads"]["create_fastq_links"]
    shell:
         " mkdir -p {output}; "
         " ln -s {input}/*.fastq.gz {output} 2>{log.std}"
