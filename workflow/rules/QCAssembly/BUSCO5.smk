#localrules: handle_busco5_output

rule busco5:
    priority: 500
    input:
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        dir=temp(directory(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{busco_lineage}.{haplotype}")),
        summary=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.summary",
        summary_json=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.summary.json",
        busco_table=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.full_table.tsv",
        missing_busco_ids=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.missing.ids",
        #summary=temp(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/short_summary.specific.{busco_lineage}.{genome_prefix}.{assembly_stage}.{haplotype}.txt"),
        #summary_json=temp(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/short_summary.specific.{busco_lineage}.{genome_prefix}.{assembly_stage}.{haplotype}.json"),
        #busco_table=temp(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/run_{busco_lineage}/full_table.tsv"),
        #missing_busco_ids=temp(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/run_{busco_lineage}/missing_busco_list.tsv"),
    #params:
    #    out_prefix= lambda wildcards: "{0}.{1}.{2}".format(wildcards.genome_prefix,
    #                                                       wildcards.assembly_stage,
    #                                                       wildcards.haplotype)
    log:
        std=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.log",
        cp=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cp.log",
        cluster_log=output_dict["cluster_log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5"],
        time=parameters["time"]["busco5"],
        mem=parameters["memory_mb"]["busco5"],
    threads:
        parameters["threads"]["busco5"]
    shell:
         " busco -m genome -l {wildcards.busco_lineage} -c {threads} -i {input.assembly} "
         " -o `basename {output.dir}` --out_path `dirname {output.dir}` > {log.std} 2>&1;"
         " cp {output.dir}/short_summary.specific.{wildcards.busco_lineage}.{wildcards.genome_prefix}.{wildcards.assembly_stage}.{wildcards.busco_lineage}.{wildcards.haplotype}.txt {output.summary} ; "
         " cp {output.dir}/short_summary.specific.{wildcards.busco_lineage}.{wildcards.genome_prefix}.{wildcards.assembly_stage}.{wildcards.busco_lineage}.{wildcards.haplotype}.json {output.summary_json} ; "
         " cp {output.dir}/run_{wildcards.busco_lineage}/full_table.tsv {output.busco_table} ; "
         " cp {output.dir}/run_{wildcards.busco_lineage}/missing_busco_list.tsv {output.missing_busco_ids} ; "


rule compress_busco5:
    priority: 500
    input:
        dir=rules.busco5.output.dir,
        #summary=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/short_summary.specific.{busco_lineage}.{genome_prefix}.{assembly_stage}.{haplotype}.txt",
        #summary_json=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/short_summary.specific.{busco_lineage}.{genome_prefix}.{assembly_stage}.{haplotype}.json",
        #busco_table=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/run_{busco_lineage}/full_table.tsv",
        #missing_busco_ids=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/run_{busco_lineage}/missing_busco_list.tsv",
    output:
        tar_gz=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.tar.gz",

    log:
        tar=output_dict["log"] / "compress_busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.tar.log",
        cluster_log=output_dict["cluster_log"] / "compress_busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "compress_busco5{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "compress_busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["compress_busco5"],
        time=parameters["time"]["compress_busco5"],
        mem=parameters["memory_mb"]["compress_busco5"],
    threads:
        parameters["threads"]["compress_busco5"]
    shell:
         " tar czf {output.tar_gz} {input.dir} > {log.tar}2>&1;"
