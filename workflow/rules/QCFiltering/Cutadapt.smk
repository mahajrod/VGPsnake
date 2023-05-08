
rule cutadapt:
    input:
        fastq=output_dict["data"] / ("fastq/{datatype}/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        fastq=output_dict["data"] / ("fastq/{datatype}/filtered/{fileprefix}%s" % config["fastq_extension"]),
        stats=output_dict["data"] / "fastq/{datatype}/filtered/{fileprefix}.cutadapt.stats"
    params:
        error_rate=lambda wildcards: "-e {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["error_rate"]) if "error_rate" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        min_read_length=lambda wildcards: " -m {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["min_read_length"]) if "min_read_length" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        min_adapter_length=lambda wildcards: " --overlap {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["min_adapter_length"]) if "min_adapter_length" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        adapter_match_times=lambda wildcards: " --times {0}".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["adapter_match_times"]) if "adapter_match_times" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        check_read_rc= lambda wildcards: " --rc " if parameters["tool_options"]["cutadapt"][wildcards.datatype]["check_read_rc"] else "",
        discard_trimmed= lambda wildcards: " --discard-trimmed " if parameters["tool_options"][wildcards.datatype]["{datatype}"]["discard_trimmed"] else "",
        forward_anywhere_adapters= lambda wildcards: " -b ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["forward_anywhere_adapter_list"]) if "forward_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        reverse_anywhere_adapters= lambda wildcards: " -B ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["reverse_anywhere_adapter_list"]) if "reverse_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        forward_three_prime_adapters= lambda wildcards: " -a ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["forward_three_prime_adapter_list"]) if "forward_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        reverse_three_prime_adapters= lambda wildcards: " -A ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["reverse_three_prime_adapter_list"]) if "reverse_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        forward_five_prime_adapters= lambda wildcards: " -g ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["forward_five_prime_adapter_list"]) if "forward_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        reverse_five_prime_adapters= lambda wildcards: " -G ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["reverse_five_prime_adapter_list"]) if "reverse_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
    log:
        std=output_dict["log"] / "cutadapt_pacbio.{datatype}.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/no_cut.cutadapt.stats.log",
        cluster_log=output_dict["cluster_log"] / "cutadapt_pacbio.{datatype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "cutadapt_pacbio.{datatype}.{fileprefix}.cluster.log"
    benchmark:
        output_dict["benchmark"] /  "cutadapt_pacbio.{datatype}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["cutadapt"],
        time=parameters["time"]["cutadapt"],
        mem=parameters["memory_mb"]["cutadapt"],
    threads:
        parameters["threads"]["cutadapt"]
    shell:
         " cutadapt -j {threads} {params.min_read_length} {params.error_rate} {params.min_adapter_length} "
         " {params.adapter_match_times} -b {params.forward_anywhere_adapters} -o {output.fastq} "
         " {params.check_read_rc} {params.discard_trimmed}"
         " {input.fastq} > {output.stats} 2>{log.std}; "

