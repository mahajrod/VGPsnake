rule merqury:
    input:
        meryl_db_dir=output_dict["kmer"] / "{0}/filtered/{0}.filtered.{1}.meryl".format(config["final_kmer_datatype"],
                                                                                        config["final_kmer_length"],),
        primary_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.hap1.fasta" % config["genome_name"]),
        alternative_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.hap2.fasta" % config["genome_name"]),
        #meryl_db_dir=(output_dict["kmer"] / "{0}/filtered/{0}.filtered.{1}.meryl".format(config["final_kmer_datatype"],
        #                                                                                config["final_kmer_length"],)).resolve(),
        #primary_assembly=(out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.hap1.p_ctg.fasta" % config["genome_name"])).resolve(),
        #alternative_assembly=(out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.hap2.p_ctg.fasta" % config["genome_name"])).resolve()
    output:
        #out_dir=output_dict["assembly_qc"] /"{assembly_stage}/merqury/{assembler}/",
        qv_file=output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.qv" % config["genome_name"]),
        completeness_stats_file=output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.completeness.stats" % config["genome_name"])
    params:
        dir=lambda wildcards: output_dict["assembly_qc"] /"{0}/merqury/{1}/".format(wildcards.assembly_stage,
                                                                                    wildcards.assembler),
        out_prefix=lambda wildcards: "{2}.{0}.{1}".format(wildcards.assembly_stage,
                                                                        wildcards.assembler,
                                                                        config["genome_name"],)
    log:
        std=output_dict["log"].resolve() / "merqury.{assembler}.{assembly_stage}.log",
        mkdir_log=(output_dict["log"]).resolve() / "merqury.{assembler}.{assembly_stage}.mkdir.log",
        cd_log=(output_dict["log"]).resolve() / "merqury.{assembler}.{assembly_stage}.cd.log",
        cluster_log=(output_dict["cluster_log"]).resolve() / "merqury.{assembler}.{assembly_stage}.cluster.log",
        cluster_err=(output_dict["cluster_error"]).resolve() / "merqury.{assembler}.{assembly_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merqury.{assembler}.{assembly_stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["merqury"],
        time=parameters["time"]["merqury"],
        mem=parameters["memory_mb"]["merqury"],
    threads:
        parameters["threads"]["merqury"]
    shell:
         " MERYL_DB=`realpath {input.meryl_db_dir}`;"
         " PRIMARY_ASSEMBLY=`realpath {input.primary_assembly}`;"
         " ALTERNATIVE_ASSEMBLY=`realpath {input.alternative_assembly}`;"
         " cd {params.dir}; "
         " OMP_NUM_THREADS={threads} merqury.sh ${{MERYL_DB}} "
         " ${{PRIMARY_ASSEMBLY}} ${{ALTERNATIVE_ASSEMBLY}} {params.out_prefix}  1>{log.std} 2>&1;"
         #" OMP_NUM_THREADS={threads} merqury.sh {input.meryl_db_dir} "
         #" {input.primary_assembly} {input.alternative_assembly} {params.out_prefix}  1>{log.std} 2>&1 || true;"


