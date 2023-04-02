rule hifiasm: # TODO: implement modes without hic data
    input:
        hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                    fileprefix=input_file_prefix_dict["hifi"],
                    allow_missing=True),
        hic_forward=input_filedict["hic"][::2],
        hic_reverse=input_filedict["hic"][1::2],
        genomescope_report=output_dict["kmer"] / ("%s/filtered/genomescope/{genome_prefix}.%s.filtered.%s.%s.genomescope.parameters" % (config["final_kmer_datatype"],
                                                                                                                                        config["final_kmer_datatype"],
                                                                                                                                        config["final_kmer_length"],
                                                                                                                                        config["final_kmer_counter"])),
    output:
        raw_unitig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.hic.r_utg.gfa",
        primary_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.hic.hap1.p_ctg.gfa",
        alternative_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.hic.hap2.p_ctg.gfa",
        primary_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hap1.gfa",
        alternative_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hap2.gfa",
        #primary_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.hifi.hic.p_ctg.gfa" % config["genome_name"]),
        #alternative_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.hifi.hic.a_ctg.gfa" % config["genome_name"]),
    params:
        purge_level=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["purge level"],
        ploidy=config["ploidy"],
        primary=lambda wildcards: "--primary" if parameters["tool_options"]["hifiasm"][wildcards.contig_options]["primary"] else "",
        output_prefix=lambda wildcards: output_dict["contig"] / "hifiasm_{0}/{1}.contig.hifi".format(wildcards.contig_options,
                                                                                                     wildcards.genome_prefix),
        cov_multiplicator=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["cov_multiplicator"],
        hic_forward=",".join(map(str, input_filedict["hic"][::2])), #in case of multiple hic libraries files in the list MUST be COMMA-separated
        hic_reverse=",".join(map(str, input_filedict["hic"][1::2]))
    log:
        std=output_dict["log"] / "hifiasm.{contig_options}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm.{contig_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`; "
         " hifiasm {params.primary} -t {threads} -l {params.purge_level}  -o {params.output_prefix} "
         " --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} "
         " --h1 {params.hic_forward} --h2 {params.hic_reverse} "
         " {input.hifi}  1>{log.std} 2>&1;"
         " ln {output.primary_contig_graph} {output.primary_alias};"
         " ln {output.alternative_contig_graph} {output.alternative_alias}"
