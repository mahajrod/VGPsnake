rule busco5:
    input:
        assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.fasta" % config["genome_name"])
    output:
        #summary=output_dict["assembly_qc"] /("{assembly_stage}/{assembler}/{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.gfa" % config["genome_name"]),
        dir=directory(output_dict["assembly_qc"] /"{assembly_stage}/busco5/{assembler}/{haplotype}"),
        summary=output_dict["assembly_qc"] /("{assembly_stage}/busco5/{assembler}/{haplotype}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}/short_summary.specific.%s.%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}.txt"
                                                   % (config["genome_name"],
                                                      config["busco_lineage"],
                                                      config["genome_name"]))
    params:
        lineage=config["busco_lineage"],
        out_prefix= lambda wildcards: "{0}.{1}.{2}.pacbio.hic.{3}".format(config["genome_name"],
                                                                          wildcards.assembly_stage,
                                                                          wildcards.assembler,
                                                                          wildcards.haplotype)
    log:
        std=output_dict["log"] / "busco5.{assembler}.{assembly_stage}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "busco5.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["busco5"],
        time=parameters["time"]["busco5"],
        mem=parameters["memory_mb"]["busco5"],
    threads:
        parameters["threads"]["busco5"]
    shell:
         " busco -m genome -l {params.lineage} -c {threads} -i {input.assembly} "
         " -o {params.out_prefix} --out_path {output.dir}  1>{log.std} 2>&1;"


#results/assembly_qc/contig/busco5/hifiasm/hap1.p/genome.contig.hifiasm.pacbio.hic.hap1.p/short_summary.specific.saccharomycetes_odb10.genome.contig.hifiasm.pacbio.hic.hap1.p.txt