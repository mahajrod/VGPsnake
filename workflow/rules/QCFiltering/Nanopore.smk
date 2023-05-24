
rule filter_nanopore:
    input:
        fastq=output_dict["data"] / ("fastq/{datatype}/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        trimmed_fastq=output_dict["data"] / ("fastq/{datatype}/trimmed/{fileprefix}%s" % config["fastq_extension"]),
        filtered_fastq=output_dict["data"] / ("fastq/{datatype}/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #stats=output_dict["data"] / "fastq/{datatype}/{stage}/{fileprefix}.cutadapt.stats"
    params:
         ab_initio=lambda wildcards: parse_option_flag("ab_initio", parameters["tool_options"]["porechop_abi"][wildcards.datatype], "--ab_initio"),
         verbosity=lambda wildcards: parse_option("verbosity", parameters["tool_options"]["porechop_abi"][wildcards.datatype], "-v"),
         porechop_abi_threads=parameters["threads"]["porechop_abi"],
         chopper_threads     =parameters["threads"]["chopper"],
         headcrop  =lambda wildcards:parse_option("headcrop",  parameters["tool_options"]["chopper"][wildcards.datatype], "--headcrop"),
         maxlength =lambda wildcards:parse_option("maxlength", parameters["tool_options"]["chopper"][wildcards.datatype], "--maxlength"),
         minlength =lambda wildcards:parse_option("minlength", parameters["tool_options"]["chopper"][wildcards.datatype], "--minlength"),
         quality   =lambda wildcards:parse_option("quality",   parameters["tool_options"]["chopper"][wildcards.datatype], "--quality"),
         tailcrop  =lambda wildcards:parse_option("tailcrop",  parameters["tool_options"]["chopper"][wildcards.datatype], "--tailcrop"),
    log:
        porechop_abi=output_dict["log"]/ "filter_nanopore.{datatype}.trimmed.{fileprefix}.porechop_abi.log",
        tee=output_dict["log"]/ "filter_nanopore.{datatype}.trimmed.{fileprefix}.tee.log",
        chopper=output_dict["log"]/ "filter_nanopore.{datatype}.trimmed.{fileprefix}.chopper.log",
        pigz_trimmed=output_dict["log"]/ "filter_nanopore.{datatype}.trimmed.{fileprefix}.pigz_trimmed.log",
        pigz_filtered=output_dict["log"]/ "filter_nanopore.{datatype}.trimmed.{fileprefix}.pigz_filtered.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "filter_nanopore.{datatype}.trimmed.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "filter_nanopore.{datatype}.trimmed.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "filter_nanopore.{datatype}.trimmed.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        cpus=parameters["threads"]["porechop_abi"] + parameters["threads"]["chopper"],
        time=parameters["time"]["filter_nanopore"],
        mem=parameters["memory_mb"]["porechop_abi"] + parameters["memory_mb"]["chopper"],
    threads:
        parameters["threads"]["porechop_abi"] + parameters["threads"]["chopper"] + 4
    shell:
        #" INPUT_FASTQ={input.fastq}; "
        #" INPUT_UNCOMPRESSED_FASTQ=${{INPUT_FASTQ%.gz}};"
        #" gunzip -c ${{INPUT_FASTQ}} > ${{}}"
        " TRIMMED_FASTQ={output.trimmed_fastq}; "
        " FILTERED_FASTQ={output.filtered_fastq}; "
        " porechop_abi {params.ab_initio} {params.verbosity} -t {params.porechop_abi_threads} "
        " -i {input.fastq} 2>{log.porechop_abi} | "
        " tee ${{TRIMMED_FASTQ%.gz}} 2>{log.tee} | "
        " chopper {params.headcrop} {params.maxlength} {params.minlength} {params.quality} {params.tailcrop} "
        " -t {params.chopper_threads} 2>{log.chopper} | pigz -p 4 > {output.filtered_fastq} 2>{log.pigz_filtered} ; "
        " pigz -p {threads} ${{TRIMMED_FASTQ%.gz}} >{log.pigz_trimmed} 2>&1; "

