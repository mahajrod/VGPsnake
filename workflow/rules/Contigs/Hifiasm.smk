rule hifiasm: # TODO: implement modes without hic data
    input:
        pacbio=expand(output_dict["data"] / ("fastq/pacbio/filtered/{fileprefix}%s" % config["fastq_extension"]),
                           fileprefix=input_file_prefix_dict["pacbio"],
                           allow_missing=True
                           ),
        hic_forward=input_filedict["hic"][::2],
        hic_reverse=input_filedict["hic"][1::2],
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{0}.filtered.{1}.{2}.genomescope.parameters".format(config["final_kmer_datatype"],
                                                                                                                               config["final_kmer_length"],
                                                                                                                               config["final_kmer_counter"]),

    output:
        raw_unitig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.pacbio.hic.r_utg.gfa" % config["genome_name"]),
        primary_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.pacbio.hic.hap1.p_ctg.gfa" % config["genome_name"]),
        alternative_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.pacbio.hic.hap2.p_ctg.gfa" % config["genome_name"]),
        #primary_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.pacbio.hic.p_ctg.gfa" % config["genome_name"]),
        #alternative_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.pacbio.hic.a_ctg.gfa" % config["genome_name"]),
    params:
        purge_level=parameters["tool_options"]["hifiasm"]["purge level"],
        ploidy=config["ploidy"],
        output_prefix=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.pacbio" % config["genome_name"]),
        cov_multiplicator=parameters["tool_options"]["hifiasm"]["cov_multiplicator"]
    log:
        std=output_dict["log"] / "hifiasm.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`; "
         " hifiasm --primary -t {threads} -l {params.purge_level}  -o {params.output_prefix} "
         " --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} "
         " --h1 {input.hic_forward} --h2 {input.hic_reverse} {input.pacbio}  1>{log.std} 2>&1;"
