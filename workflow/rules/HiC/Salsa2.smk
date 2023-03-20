
rule salsa2: #
    input:
        bed=out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/alignment/%s.hic_scaffolding.{assembler}.{haplotype}.bwa.filtered.rmdup.bed"  % config["genome_name"]),
        reference=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"]),
        reference_fai=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fai" % config["genome_name"])
    output:
        dir=directory(out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/scaffolding/%s"  % config["genome_name"])),
    params:
        min_contig_len=parameters["tool_options"]["salsa2"]["min_contig_len"],
        restriction_seq=parameters["tool_options"]["salsa2"]["restriction_seq"][config["hic_enzyme_set"]],
        sortorder=parameters["tool_options"]["pretextmap"]["sortorder"]
    log:
        std=output_dict["log"]  / "salsa2.{assembler}.hic_scaffolding.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "salsa2.{assembler}.hic_scaffolding.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "salsa2.{assembler}.hic_scaffolding.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "salsa2.{assembler}.hic_scaffolding.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["salsa2"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["salsa2"]["yaml"])
    resources:
        cpus=parameters["threads"]["salsa2"] ,
        time=parameters["time"]["salsa2"],
        mem=parameters["memory_mb"]["salsa2"]
    threads: parameters["threads"]["salsa2"]

    shell:
        " run_pipeline.py -a {input.reference} -l {input.reference_fai} -b {input.bed} -e {params.restriction_seq} "
        " -m yes -p yes "
        " -o {output.dir}  > {log.std} 2>&1"
