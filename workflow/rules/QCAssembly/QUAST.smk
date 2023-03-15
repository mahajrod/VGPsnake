rule quast:
    input:
        #primary_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.p_ctg.fasta" % config["genome_name"]),
        #alternative_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.a_ctg.fasta" % config["genome_name"])
        assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.fasta" % config["genome_name"]),
    output:
        #summary=output_dict["assembly_qc"] /("{assembly_stage}/{assembler}/{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.gfa" % config["genome_name"]),
        dir=directory(output_dict["assembly_qc"] /("{assembly_stage}/quast/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}"
                                                   % config["genome_name"]))
    params:
        large_genome_flag="--large" if parameters["tool_options"]["quast"]["large_genome"] else "",
    log:
        std=output_dict["log"] / "quast.{assembler}.{assembly_stage}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "quast.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "quast.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "quast.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["quast"],
        time=parameters["time"]["quast"],
        mem=parameters["memory_mb"]["quast"],
    threads:
        parameters["threads"]["quast"]
    shell:
         " quast -o {output.dir} -t {threads} -l {wildcards.haplotype} {params.large_genome_flag} "
         " {input.assembly}  1>{log.std} 2>&1;"
