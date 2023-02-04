
rule fastqc:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        fastq=output_dict["data"] / ("fastq/{datatype}/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        zip=output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fastqc.stats"
    params:
        kmer=parameters["tool_options"]["fastqc"]["kmer_length"],
        out_dir=lambda wildcards: output_dict["qc"] / "fastqc/{0}/{1}/".format(wildcards.datatype,
                                                                               wildcards.stage),
        nogroup=lambda wildcards: "" if wildcards.datatype in config["long_read_data"] else "--nogroup" # turns off base grouping for short reads

    log:
        std=output_dict["log"]/ "fastqc.{datatype}.{stage}.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"]/ "fastqc.{datatype}.{stage}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "fastqc.{datatype}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "fastqc.{datatype}.{stage}.{fileprefix}.benchmark.txt"
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
        " fastqc {params.nogroup} -k {params.kmer} -t {threads} -o {params.out_dir} {input} 1>{log.std} 2>&1; "
        #" workflow/scripts/convert_fastqc_output.py -f {output.forward_fastqc} -r {output.reverse_fastqc} "
        #" -s {wildcards.library_id} -o {output.stats} 1>{log.stats} 2>&1 "

rule multiqc:
    input:
        output_dict["qc"] / "fastqc/{datatype}/{stage}/"
    output:
        dir=directory(output_dict["qc"] / "multiqc/{datatype}/{stage}/"),
        report=output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report"
        #stats=merged_raw_multiqc_dir_path / "{library_id}/{library_id}.raw.multiqc.stats"

    log:
        std=output_dict["log"]/ "multiqc.{datatype}.{stage}.log",
        #stats=log_dir_path / "{library_id}/multiqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"]/ "multiqc.{datatype}.{stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "multiqc.{datatype}.{stage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "multiqc.{datatype}.{stage}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["multiqc"],
        time=parameters["time"]["multiqc"],
        mem=parameters["memory_mb"]["multiqc"],
    threads:
        parameters["threads"]["multiqc"]
    shell:
        " multiqc --filename {output.report} -p --outdir {output.dir} --comment {wildcards.datatype} {input}  2>&1; "