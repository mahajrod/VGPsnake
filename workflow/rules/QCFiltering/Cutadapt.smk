
rule cutadapt:
    input:
        fastq=output_dict["data"] / ("fastq/hifi/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        fastq=output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        stats=output_dict["data"] / "fastq/hifi/filtered/{fileprefix}.cutadapt.stats"
    params:
        error_rate=lambda wildcards: "-e {0} ".format(parameters["tool_options"]["cutadapt"]["hifi"]["error_rate"]) if "error_rate" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        min_read_length=lambda wildcards: " -m {0} ".format(parameters["tool_options"]["cutadapt"]["hifi"]["min_read_length"]) if "min_read_length" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        min_adapter_length=lambda wildcards: " --overlap {0} ".format(parameters["tool_options"]["cutadapt"]["hifi"]["min_adapter_length"]) if "min_adapter_length" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        adapter_match_times=lambda wildcards: " --times {0}".format(parameters["tool_options"]["cutadapt"]["hifi"]["adapter_match_times"]) if "adapter_match_times" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        check_read_rc= lambda wildcards: " --rc " if ( ("check_read_rc" in parameters["tool_options"]["cutadapt"]["hifi"]) and  parameters["tool_options"]["cutadapt"]["hifi"]["check_read_rc"]) else "",
        discard_trimmed= lambda wildcards: " --discard-trimmed " if ( ("discard_trimmed" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["discard_trimmed"]) else "",
        forward_anywhere_adapters= lambda wildcards: " -b ".join(parameters["tool_options"]["cutadapt"]["hifi"]["forward_anywhere_adapter_list"]) if "forward_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        reverse_anywhere_adapters= lambda wildcards: " -B ".join(parameters["tool_options"]["cutadapt"]["hifi"]["reverse_anywhere_adapter_list"]) if "reverse_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        forward_three_prime_adapters= lambda wildcards: " -a ".join(parameters["tool_options"]["cutadapt"]["hifi"]["forward_three_prime_adapter_list"]) if "forward_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        reverse_three_prime_adapters= lambda wildcards: " -A ".join(parameters["tool_options"]["cutadapt"]["hifi"]["reverse_three_prime_adapter_list"]) if "reverse_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        forward_five_prime_adapters= lambda wildcards: " -g ".join(parameters["tool_options"]["cutadapt"]["hifi"]["forward_five_prime_adapter_list"]) if "forward_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        reverse_five_prime_adapters= lambda wildcards: " -G ".join(parameters["tool_options"]["cutadapt"]["hifi"]["reverse_five_prime_adapter_list"]) if "reverse_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
    log:
        std=output_dict["log"] / "cutadapt_pacbio.hifi.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/no_cut.cutadapt.stats.log",
        cluster_log=output_dict["cluster_log"] / "cutadapt_pacbio.hifi.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "cutadapt_pacbio.hifi.{fileprefix}.cluster.log"
    benchmark:
        output_dict["benchmark"] /  "cutadapt_pacbio.hifi.{fileprefix}.benchmark.txt"
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
         " {params.adapter_match_times} {params.forward_anywhere_adapters} -o {output.fastq} "
         " {params.check_read_rc} {params.discard_trimmed} "
         " {input.fastq} > {output.stats} 2>{log.std}; "


rule cutadapt_paired:
    input:
        forward_fastq=lambda wildcards: output_dict["data"] / ("fastq/{0}/raw/{1}{2}{3}".format("illumina",
                                                                                                wildcards.pairprefix,
                                                                                                input_forward_suffix_dict["illumina"],
                                                                                                config["fastq_extension"])),
        reverse_fastq=lambda wildcards: output_dict["data"] / ("fastq/{0}/raw/{1}{2}{3}".format("illumina",
                                                                                                wildcards.pairprefix,
                                                                                                input_reverse_suffix_dict["illumina"],
                                                                                                config["fastq_extension"])),
        #lambda wildcards: input_filedict["illumina"][::2] if "illumina" in config["paired_fastq_based_data"] else [],
        #reverse_fastq=lambda wildcards: input_filedict["illumina"][1::2] if "illumina" in config["paired_fastq_based_data"] else [],
        #fastq=output_dict["data"] / ("fastq/illumina/raw/{pairprefix}%s" % config["fastq_extension"])
    output:
        #fastq=output_dict["data"] / ("fastq/illumina/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #forward_fastq=output_dict["data"] / ("fastq/illumina/filtered/{pairprefix}.paired_1%s" % config["fastq_extension"]),
        forward_fastq=output_dict["data"] / ("fastq/illumina/filtered/{pairprefix}_1%s" % config["fastq_extension"]),
        reverse_fastq=output_dict["data"] / ("fastq/illumina/filtered/{pairprefix}_2%s" % config["fastq_extension"]),
        #reverse_fastq=output_dict["data"] / ("fastq/illumina/filtered/{pairprefix}.paired_2%s" % config["fastq_extension"]),
        stats=output_dict["data"] / "fastq/illumina/filtered/{pairprefix}.cutadapt.stats"
    params:
        error_rate=lambda wildcards: "-e {0} ".format(parameters["tool_options"]["cutadapt"]["illumina"]["error_rate"]) if "error_rate" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        min_read_length=lambda wildcards: " -m {0} ".format(parameters["tool_options"]["cutadapt"]["illumina"]["min_read_length"]) if "min_read_length" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        min_adapter_length=lambda wildcards: " --overlap {0} ".format(parameters["tool_options"]["cutadapt"]["illumina"]["min_adapter_length"]) if "min_adapter_length" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        adapter_match_times=lambda wildcards: " --times {0}".format(parameters["tool_options"]["cutadapt"]["illumina"]["adapter_match_times"]) if "adapter_match_times" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        check_read_rc= lambda wildcards: " --rc " if ( ("check_read_rc" in parameters["tool_options"]["cutadapt"]["illumina"]) and  parameters["tool_options"]["cutadapt"]["illumina"]["check_read_rc"]) else "",
        discard_trimmed= lambda wildcards: " --discard-trimmed " if ( ("discard_trimmed" in parameters["tool_options"]["cutadapt"]["illumina"]) and parameters["tool_options"]["cutadapt"]["illumina"]["discard_trimmed"]) else "",
        forward_anywhere_adapters= lambda wildcards: " -b ".join(parameters["tool_options"]["cutadapt"]["illumina"]["forward_anywhere_adapter_list"]) if "forward_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        reverse_anywhere_adapters= lambda wildcards: " -B ".join(parameters["tool_options"]["cutadapt"]["illumina"]["reverse_anywhere_adapter_list"]) if "reverse_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        forward_three_prime_adapters= lambda wildcards: " -a ".join(parameters["tool_options"]["cutadapt"]["illumina"]["forward_three_prime_adapter_list"]) if "forward_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        reverse_three_prime_adapters= lambda wildcards: " -A ".join(parameters["tool_options"]["cutadapt"]["illumina"]["reverse_three_prime_adapter_list"]) if "reverse_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        forward_five_prime_adapters= lambda wildcards: " -g ".join(parameters["tool_options"]["cutadapt"]["illumina"]["forward_five_prime_adapter_list"]) if "forward_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
        reverse_five_prime_adapters= lambda wildcards: " -G ".join(parameters["tool_options"]["cutadapt"]["illumina"]["reverse_five_prime_adapter_list"]) if "reverse_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["illumina"] else "",
    log:
        std=output_dict["log"] / "cutadapt_pacbio.illumina.{pairprefix}.log",
        #stats=log_dir_path / "{library_id}/no_cut.cutadapt.stats.log",
        cluster_log=output_dict["cluster_log"] / "cutadapt_pacbio.illumina.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "cutadapt_pacbio.illumina.{pairprefix}.cluster.log"
    benchmark:
        output_dict["benchmark"] /  "cutadapt_pacbio.illumina.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["cutadapt"],
        time=parameters["time"]["cutadapt"],
        mem=parameters["memory_mb"]["cutadapt"],
    threads:
        parameters["threads"]["cutadapt"]
    shell:
         " cutadapt --paired -j {threads} {params.min_read_length} {params.error_rate} {params.min_adapter_length} "
         " {params.adapter_match_times} {params.forward_anywhere_adapters} "
         " {params.check_read_rc} {params.discard_trimmed} "
         " -o {output.forward_fastq} -p {output.reverse_fastq} "
         " {input.forward_fastq} {input.reverse_fastq} > {output.stats} 2>{log.std}; "


