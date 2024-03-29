"""
rule minimap2_index:
    input:
        reference=output_dict["contig"] / ("{assembler}/%s.purge_dups.{assembler}.hifi.hic.{haplotype}_ctg.fasta" % config["genome_name"])
    priority: 1000
    output:
        index=output_dict["contig"] / ("{assembler}/%s.purge_dups.{assembler}.hifi.hic.{haplotype}_ctg.minimap2.idx" % config["genome_name"])
    params:
        index_size=parameters["tool_options"]["minimap2"]["index_size"],
    log:
        minimap2_index=output_dict["log"] / "minimap2_index.{assembler}.purge_dups.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_index.{assembler}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_index.{assembler}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_index.{assembler}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ../../../config["conda"]["common"]["yaml"]
    resources:
        cpus=parameters["threads"]["minimap2_index"],
        time=parameters["time"]["minimap2_index"],
        mem=parameters["memory_mb"]["minimap2_index"]
    threads: parameters["threads"]["minimap2_index"]
    shell:
        " minimap2 -t {threads} -I {params.index_size} -d {output.index} {input.reference} > {log.minimap2_index} 2>&1"
"""
localrules: create_primary_contig_link, create_link_for_purged_fasta, merge_pri_hapdups_with_alt
ruleorder: create_primary_contig_link > merge_pri_hapdups_with_alt

rule create_primary_contig_link:
    input:
        #fasta=out_dir_path / ("contig/{assembler}/%s.contig.{assembler}.hap1.fasta" % config["genome_name"]),
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.hap1.fasta" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                            stage_dict["purge_dups"]["prev_stage"]))
    output:
        fasta=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.hap1.fasta"
        #fasta=out_dir_path / ("purge_dups/{assembler}/input/%s.contig.{assembler}.hap1.fasta" % config["genome_name"])
    log:
        std=output_dict["log"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.log",
        cluster_log=output_dict["cluster_log"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln {input.fasta} {output.fasta} 1>{log.std} 2>&1"

rule minimap2_purge_dups_reads: # TODO: add nanopore support
    input:
        fastq=output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        reference=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.fasta"
    output:
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.{haplotype}.{fileprefix}.paf.gz"
        #paf=out_dir_path  / ("purge_dups/{assembler}/{haplotype}/%s.purge_dups.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        index_size=parameters["tool_options"]["minimap2"]["index_size"],
        mapping_scheme=parameters["tool_options"]["minimap2"]["hifi_alignment_scheme"], # TODO: make this adjustable depending on read type
    log:
        std=output_dict["log"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " minimap2 {params.mapping_scheme} -I {params.index_size} -t {threads}  {input.reference} "
        " {input.fastq} 2>{log.std} |  gzip -c - > {output.paf} "

rule get_purge_dups_read_stat: #TODO: adjust -d -m -u options for calcuts
    input:
        #paf=expand(out_dir_path / ("purge_dups/{assembler}/{haplotype}/%s.purge_dups.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"]),
        #                   fileprefix=input_file_prefix_dict["hifi"],
        #                   allow_missing=True),
        paf=expand(out_dir_path / ("purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/%s.{haplotype}.{fileprefix}.paf.gz" % config["genome_prefix"]),
                           fileprefix=input_file_prefix_dict["hifi"],
                           allow_missing=True),
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{1}.{0}.filtered.{2}.{3}.genomescope.parameters".format(config["final_kmer_datatype"],
                                                                                                                                   config["genome_prefix"],
                                                                                                                                   config["final_kmer_length"],
                                                                                                                                   config["final_kmer_counter"])
    output:
        pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.stat",
        pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.base.cov",
        cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/cutoffs"
    params:
        out_dir=lambda wildcards: out_dir_path  / "purge_dups/{0}..{1}/{2}".format(wildcards.prev_stage_parameters,
                                                                                  wildcards.purge_dups_parameters,
                                                                                  wildcards.haplotype),
        cov_multiplicator=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["cov_multiplicator"]
        #cov_multiplicator=parameters["tool_options"]["purge_dups"]["cov_multiplicator"]

    log:
        pbstat=output_dict["log"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.pbstat.log",
        calcuts=output_dict["log"]  / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.calcuts.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["get_purge_dups_read_stat"] ,
        time=parameters["time"]["get_purge_dups_read_stat"],
        mem=parameters["memory_mb"]["get_purge_dups_read_stat"]
    threads: parameters["threads"]["get_purge_dups_read_stat"]

    shell:
        " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`;"
        " pbcstat -O {params.out_dir} {input.paf} 1>{log.pbstat} 2>&1; "
        " calcuts -d 1 -u ${{COV_UPPER_BOUNDARY}} {output.pbstat} > {output.cutoffs} 2>{log.calcuts} " #check parameters for calcuts

rule minimap2_purge_dups_assembly:
    input:
        reference = out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.fasta"
        #reference=out_dir_path  / ("purge_dups/{assembler}/input/%s.contig.{assembler}.{haplotype}.fasta" % config["genome_name"])
    output:
        split_reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.purge_dups_input.{haplotype}.split.fasta",
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.purge_dups_input.{haplotype}.split.minimap2.self.paf.gz"
        #paf=output_dict["purge_dups"] / ("{assembler}/%s.purge_dups.{assembler}.hifi.hic.{haplotype}_ctg.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        index_size=parameters["tool_options"]["minimap2"]["index_size"],
        mapping_scheme=parameters["tool_options"]["minimap2"]["self_alignment_scheme"]
    log:
        split_fa=output_dict["log"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.split_fa.log",
        minimap2=output_dict["log"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.minimap2.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " split_fa {input.reference} > {output.split_reference} 2>{log.split_fa};"
        " minimap2 {params.mapping_scheme} -I {params.index_size} -t {threads}  {output.split_reference} "
        " {output.split_reference} 2>{log.minimap2} |  gzip -c - > {output.paf} "

rule purge_dups: # TODO: find what options are used in ERGA for get_seqs
    input:
        cutoffs=rules.get_purge_dups_read_stat.output.cutoffs,
        pbbasecov=rules.get_purge_dups_read_stat.output.pbbasecov,
        self_paf=rules.minimap2_purge_dups_assembly.output.paf,
        reference = out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.fasta"
    output:
        bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.dups.bed",
        purged=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.purge_dups.{haplotype}.purged.fasta",
        hapdups=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.purge_dups.{haplotype}.hap.fasta",
    params:
        #bed_local_path=lambda wildcards: "{0}.dups.bed".format(wildcards.genome_prefix),
        out_dir=lambda wildcards: out_dir_path  / "purge_dups/{0}..{1}/{2}".format(wildcards.prev_stage_parameters,
                                                                                 wildcards.purge_dups_parameters,
                                                                                 wildcards.haplotype),
        get_seq_prefix=lambda wildcards: "{0}.purge_dups.{1}".format(wildcards.genome_prefix, wildcards.haplotype)
    log:
        purge_dups=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.purge_dups.log",
        get_seqs=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.get_seqs.log",
        cluster_log=output_dict["cluster_log"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]

    shell:
        " purge_dups -2 -T {input.cutoffs} -c {input.pbbasecov} {input.self_paf} > {output.bed} 2>{log.purge_dups};"
        " PURGE_DUPS_BED=`realpath {output.bed}`;"
        " REFERENCE=`realpath {input.reference}`;"
        " GET_SEQ_LOG=`realpath {log.get_seqs}`;"
        " cd {params.out_dir};"
        " get_seqs -p {params.get_seq_prefix} ${{PURGE_DUPS_BED}} ${{REFERENCE}} > ${{GET_SEQ_LOG}} 2>&1;"
        " for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done"
        #" cp dups.bed {params.bed_local_path} "


rule merge_pri_hapdups_with_alt: # TODO: add handling of polyploid cases
    input:
        alt_contig=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.hap2.fasta" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                                 stage_dict["purge_dups"]["prev_stage"])),
        #reference=out_dir_path  / ("contig/{assembler}/%s.contig.{assembler}.hap2.fasta" % config["genome_name"]),
        pri_hapdups=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/hap1/{genome_prefix}.purge_dups.hap1.hap.fasta",
        #pri_hapdups=out_dir_path / ("purge_dups/{assembler}/hap1/%s.purge_dups.{assembler}.hap1.hap.fasta" % config["genome_name"])
    output:
        alt_plus_pri_hapdup=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.hap2.fasta",
        #alt_plus_pri_hapdup=out_dir_path  / ("purge_dups/{assembler}/input/%s.contig.{assembler}.hap2.fasta" % config["genome_name"]),
    log:
        std=output_dict["log"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.log",
        cluster_log=output_dict["cluster_log"] / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["merge_pri_hapdups_with_alt"] ,
        time=parameters["time"]["merge_pri_hapdups_with_alt"],
        mem=parameters["memory_mb"]["merge_pri_hapdups_with_alt"]
    threads: parameters["threads"]["merge_pri_hapdups_with_alt"]

    shell:
        " cat {input.alt_contig} {input.pri_hapdups} > {output.alt_plus_pri_hapdup} 2>{log.std}"

rule create_link_for_purged_fasta:
    input:
        purged=rules.purge_dups.output.purged
    output:
        purged=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta"
    log:
        std=output_dict["log"]  / "create_link_for_purged_fasta.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "create_link_for_purged_fasta.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_link_for_purged_fasta.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln {input.purged} {output.purged} > {log.std} 2>&1;"

