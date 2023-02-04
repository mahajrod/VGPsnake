localrules: create_fastq_links

rule create_fastq_links:
    priority: 1000
    input:
        input_dir_path.resolve() / "{datatype}/fastq/{fileprefix}" %  config["fastq_extension"]
    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        output_dict["data"] / "fastq/{datatype}/raw/{fileprefix}%s" % config["fastq_extension"]
    log:
        std=output_dict["log"] / "preprocessing/{datatype}/create_fastq_links.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "preprocessing/{datatype}/create_fastq_links.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "preprocessing/{datatype}/create_fastq_links.{fileprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "preprocessing/{datatype}/create_fastq_links.{fileprefix}.benchmark.txt",
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
         " ln -s {input} {output} 2>{log.std}"
