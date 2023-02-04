
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
        " multiqc --filename {output.report} -p --outdir {output.dir} --comment {wildcards.datatype} {input}>{log.std} 2>&1; "