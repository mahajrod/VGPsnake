rule busco5:
    input:
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype, [^.]+}.fasta"
    output:
        dir=directory(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype, [^.]+}"),
    params:
        lineage=config["busco_lineage"],
        out_prefix= lambda wildcards: "{0}.{1}.{2}".format(wildcards.genome_prefix,
                                                           wildcards.assembly_stage,
                                                           wildcards.haplotype)
    log:
        std=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype, [^.]+}.log",
        cluster_log=output_dict["cluster_log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype, [^.]+}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype, [^.]+}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype, [^.]+}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5"],
        time=parameters["time"]["busco5"],
        mem=parameters["memory_mb"]["busco5"],
    threads:
        parameters["threads"]["busco5"]
    shell:
         " busco -m genome -l {params.lineage} -c {threads} -i {input.assembly} "
         " -o {params.out_prefix} --out_path {output.dir}  1>{log.std} 2>&1;"

