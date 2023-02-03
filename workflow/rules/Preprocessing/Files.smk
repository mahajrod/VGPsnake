localrules: create_fastq_links

rule create_fastq_links:
    priority: 1000
    input:
        directory(input_dir_path / "{datatype}/fastq")

    output:
        directory(output_dict["data"] / "{datatype}/fastq/")
    log:
        std=log_dir_path / "preprocessing/{datatype}/create_fastq_links.log",
        cluster_log=cluster_log_dir_path / "preprocessing/{datatype}/create_fastq_links.cluster.log",
        cluster_err=cluster_log_dir_path / "preprocessing/{datatype}/create_fastq_links.cluster.err",
    benchmark:
        benchmark_dir_path / "preprocessing/{datatype}/create_fastq_links.benchmark.txt",
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
