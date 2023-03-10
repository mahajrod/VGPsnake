rule merqury:
    input:
        meryl_db_dir=output_dict["kmer"] / "{0}/filtered/{0}.filtered.{1}.meryl".format(config["final_kmer_datatype"],
                                                                                        config["final_kmer_length"],) ,
        primary_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.p_ctg.fasta" % config["genome_name"]),
        alternative_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.a_ctg.fasta" % config["genome_name"])
    output:
        #summary=output_dict["assembly_qc"] /("{assembly_stage}/{assembler}/{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.gfa" % config["genome_name"]),
        dir=directory(output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic"
                                                   % config["genome_name"]))
    params:
        out_prefix=lambda wildcards: output_dict["assembly_qc"] /("{0}/merqury/{1}/{2}.{0}.{1}.pacbio.hic".format(wildcards.assembly_stage,
                                                                                                                  wildcards.assembler,
                                                                                                                  config["genome_name"],))
    log:
        std=output_dict["log"] / "merqury.{assembler}.{assembly_stage}.log",
        cluster_log=output_dict["cluster_log"] / "merqury.{assembler}.{assembly_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merqury.{assembler}.{assembly_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merqury.{assembler}.{assembly_stage}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["merqury"],
        time=parameters["time"]["merqury"],
        mem=parameters["memory_mb"]["merqury"],
    threads:
        parameters["threads"]["merqury"]
    shell:
         " OMP_NUM_THREADS={threads} merqury.sh {input.meryl_db_dir} "
         " {input.primary_assembly} {input.alternative_assembly} {params.out_prefix}  1>{log.std} 2>&1;"

