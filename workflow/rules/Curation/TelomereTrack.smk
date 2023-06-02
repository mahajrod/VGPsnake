localrules: create_curation_input_links

rule telo_finder: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta",
    output:
        res=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.canonical.txt",
    params:
        size=parse_option("size", parameters["tool_options"]["telo_finder"],  "--size", default_value="default"),
        min_kmer=parse_option("min_kmer", parameters["tool_options"]["telo_finder"], "--klo", default_value="default"),
        max_kmer=parse_option("max_kmer", parameters["tool_options"]["telo_finder"], "--khi", default_value="default"),
        ends=parse_option("ends", parameters["tool_options"]["telo_finder"], "--ends", default_value="default")
    log:
        std=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.log",
        ln=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cp=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cp.log",
        cluster_log=output_dict["cluster_log"] / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"]
    threads: parameters["threads"]["telo_finder"]

    shell:
        " WORKDIR=`dirname {output.res}`; "
        " cd ${{WORKDIR}}; "
        " workflow/script/rapid_curation/telo_finder.py {params.size} {params.max_kmer} {params.max_kmer} "
        " {params.ends} {input.fasta} > {log.std} 2>&1; "
        " cp canonical.txt `basename {output.res}` > {log.cp} 2>&1; "