
rule fastqc:
    input:
        fastq=output_dict["data"]  / "{datatype}/fastq/{fileprefix}.fastq.gz"
    output:
        #dir=directory(merged_raw_fastqc_dir_path / "{library_id}"),
        zip=output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fastqc.stats"
    params:
        kmer=parameters["fastqc"]["kmer_length"],
        out_dir=lambda wildcards: output_dict["qc"] / "fastqc/{0}/{1}/".format(wildcards.datatype,
                                                                               wildcards.stage),
    log:
        std=output_dict["log"]/ "fastqc/{datatype}/{stage}/fastqc.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"]/ "fastqc/{datatype}/{stage}/fastqc.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "fastqc/{datatype}/{stage}/fastqc.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "fastqc/{datatype}/{stage}/fastqc.{fileprefix}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["fastqc"],
        time=parameters["time"]["fastqc"],
        mem=parameters["memory_mb"]["fastqc"],
    threads:
        parameters["threads"]["fastqc"]
    shell:
        #" mkdir -p {output.dir}; "
        " fastqc --nogroup -k {params.kmer} -t {threads} -o {params.out_dir} {input} 1>{log.std} 2>&1; "
        #" workflow/scripts/convert_fastqc_output.py -f {output.forward_fastqc} -r {output.reverse_fastqc} "
        #" -s {wildcards.library_id} -o {output.stats} 1>{log.stats} 2>&1 "

