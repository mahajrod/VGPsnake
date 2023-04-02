
rule cutadapt_pacbio:
    input:
        fastq=output_dict["data"] / ("fastq/hifi/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        fastq=output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        stats=output_dict["data"] / "fastq/hifi/filtered/{fileprefix}.cutadapt.stats"

    params:
        error_rate=parameters["tool_options"]["cutadapt"]["hifi"]["error_rate"],
        min_adapter_length=parameters["tool_options"]["cutadapt"]["hifi"]["min_adapter_length"],
        adapter_match_times=parameters["tool_options"]["cutadapt"]["hifi"]["adapter_match_times"],
        check_read_rc= " --rc" if parameters["tool_options"]["cutadapt"]["hifi"]["check_read_rc"] else "",
        discard_trimmed= " --discard-trimmed" if parameters["tool_options"]["cutadapt"]["hifi"]["discard_trimmed"] else "",
        anywhere_adapters= " -b ".join(parameters["tool_options"]["cutadapt"]["hifi"]["anywhere_adapter_list"])
    log:
        std=output_dict["log"] / "cutadapt_pacbio.hifi.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/no_cut.cutadapt.stats.log",
        cluster_log=output_dict["cluster_log"] / "cutadapt_pacbio.hifi.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "cutadapt_pacbio.hifi.{fileprefix}.cluster.log"
    benchmark:
        output_dict["benchmark"] /  "cutadapt_pacbio.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["cutadapt"],
        time=parameters["time"]["cutadapt"],
        mem=parameters["memory_mb"]["cutadapt"],
    threads:
        parameters["threads"]["cutadapt"]
    shell:
         " cutadapt -j {threads} -e {params.error_rate} --overlap {params.min_adapter_length} "
         " --times {params.adapter_match_times} -b {params.anywhere_adapters} -o {output.fastq} "
         " {params.check_read_rc} {params.discard_trimmed}"
         " {input.fastq} > {output.stats} 2>{log.std}; "

