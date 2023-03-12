"""
rule minimap2_index:
    input:
        reference=output_dict["contig"] / ("{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.fasta" % config["genome_name"])
    priority: 1000
    output:
        index=output_dict["contig"] / ("{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.minimap2.idx" % config["genome_name"])
    params:
        index_size=parameters["tool_options"]["minimap2"]["index_size"],
    log:
        minimap2_index=output_dict["log"] / "minimap2_index.{assembler}.{assembly_stage}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_index.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_index.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_index.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["minimap2_index"],
        time=parameters["time"]["minimap2_index"],
        mem=parameters["memory_mb"]["minimap2_index"]
    threads: parameters["threads"]["minimap2_index"]
    shell:
        " minimap2 -t {threads} -I {params.index_size} -d {output.index} {input.reference} > {log.minimap2_index} 2>&1"
"""
localrules: create_primary_contig_link, create_link_for_purged_fasta
ruleorder: create_primary_contig_link > merge_pri_hapdups_with_alt

rule create_primary_contig_link:
    input:
        fasta=out_dir_path / ("contig/{assembler}/%s.contig.{assembler}.hap1.fasta" % config["genome_name"]),
    output:
        fasta=out_dir_path / ("{assembly_stage}/{assembler}/input/%s.contig.{assembler}.hap1.fasta" % config["genome_name"])
    log:
        std=output_dict["log"]  / "create_contig_links.{assembly_stage}.{assembler}.hap1.log",
        cluster_log=output_dict["cluster_log"] / "create_contig_links.{assembly_stage}.{assembler}.hap1.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_contig_links.{assembly_stage}.{assembler}.hap1.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_contig_links.{assembly_stage}.{assembler}.hap1.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln {input.fasta} {output.fasta} 1>{log.std} 2>&1"

rule minimap2_purge_dups_reads: # TODO: add nanopore support
    input:
        fastq=output_dict["data"] / ("fastq/pacbio/filtered/{fileprefix}%s" % config["fastq_extension"]),
        reference=out_dir_path  / ("{assembly_stage}/{assembler}/input/%s.contig.{assembler}.{haplotype}.fasta" % config["genome_name"])
    output:
        paf=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/%s.{assembly_stage}.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        index_size=parameters["tool_options"]["minimap2"]["index_size"],
        mapping_scheme=parameters["tool_options"]["minimap2"]["hifi_alignment_scheme"], # TODO: make this adjustable depending on read type
    log:
        std=output_dict["log"]  / "minimap2_purge_dups_reads.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_reads.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_reads.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{assembler}.{assembly_stage}.{haplotype}.{fileprefix}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
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
        paf=expand(out_dir_path / ("{assembly_stage}/{assembler}/{haplotype}/%s.{assembly_stage}.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"]),
                           fileprefix=input_file_prefix_dict["pacbio"],
                           allow_missing=True),
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{0}.filtered.{1}.{2}.genomescope.parameters".format(config["final_kmer_datatype"],
                                                                                                                               config["final_kmer_length"],
                                                                                                                               config["final_kmer_counter"])
    output:
        pbstat=out_dir_path /  "{assembly_stage}/{assembler}/{haplotype}/PB.stat",
        pbbasecov=out_dir_path /  "{assembly_stage}/{assembler}/{haplotype}/PB.base.cov",
        cutoffs=out_dir_path /  "{assembly_stage}/{assembler}/{haplotype}/cutoffs"
    params:
        out_dir=lambda wildcards: out_dir_path  / "{0}/{1}/{2}".format(wildcards.assembly_stage, wildcards.assembler, wildcards.haplotype),
        cov_multiplicator=parameters["tool_options"]["hifiasm"]["cov_multiplicator"]

    log:
        pbstat=output_dict["log"] / "get_purge_dups_read_stat.{assembler}.{assembly_stage}.{haplotype}.pbstat.log",
        calcuts=output_dict["log"]  / "get_purge_dups_read_stat.{assembler}.{assembly_stage}.{haplotype}.calcuts.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_dups_read_stat.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_dups_read_stat.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
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
        reference=out_dir_path  / ("{assembly_stage}/{assembler}/input/%s.contig.{assembler}.{haplotype}.fasta" % config["genome_name"])
    output:
        split_reference=out_dir_path / ("{assembly_stage}/{assembler}/{haplotype}/%s.{assembly_stage}.{assembler}.{haplotype}.split.fasta" % config["genome_name"]),
        paf=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/%s.{assembly_stage}.{assembler}.{haplotype}.split.minimap2.self.paf.gz" % config["genome_name"])
        #paf=output_dict["purge_dups"] / ("{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        index_size=parameters["tool_options"]["minimap2"]["index_size"],
        mapping_scheme=parameters["tool_options"]["minimap2"]["self_alignment_scheme"]
    log:
        split_fa=output_dict["log"]  / "minimap2_purge_dups_assembly.{assembler}.{assembly_stage}.{haplotype}.split_fa.log",
        minimap2=output_dict["log"]  / "minimap2_purge_dups_assembly.{assembler}.{assembly_stage}.{haplotype}.minimap2.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_assembly.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_assembly.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_assembly.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
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
        reference=out_dir_path  / ("{assembly_stage}/{assembler}/input/%s.contig.{assembler}.{haplotype}.fasta" % config["genome_name"])
    output:
        bed=out_dir_path  / "{assembly_stage}/{assembler}/{haplotype}/dups.bed",
        purged=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/%s.{assembly_stage}.{assembler}.{haplotype}.purged.fasta" % config["genome_name"]),
        hapdups=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/%s.{assembly_stage}.{assembler}.{haplotype}.hap.fasta" % config["genome_name"]),
    params:
        out_dir=lambda wildcards: out_dir_path  / "{0}/{1}/{2}".format(wildcards.assembly_stage,
                                                                       wildcards.assembler,
                                                                       wildcards.haplotype),
        get_seq_prefix=lambda wildcards: "{1}.{2}.{0}.{3}".format(wildcards.assembler,
                                                                  config["genome_name"],
                                                                  wildcards.assembly_stage,
                                                                   wildcards.haplotype)
    log:
        purge_dups=output_dict["log"]  / "purge_dups.{assembler}.{assembly_stage}.{haplotype}.purge_dups.log",
        get_seqs=output_dict["log"]  / "purge_dups.{assembler}.{assembly_stage}.{haplotype}.get_seqs.log",
        cluster_log=output_dict["cluster_log"] / "purge_dups.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "purge_dups.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "purge_dups.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]

    shell:
        " purge_dups -2 -T {input.cutoffs} -c {input.pbbasecov} {input.self_paf} > {output.bed} 2>{log.purge_dups};"
        " PURGE_DUPS_BED=`realpath {output.bed}`;"
        " REFERENCE=`realpath {input.reference}`;"
        " cd {params.out_dir};"
        " get_seqs -p {params.get_seq_prefix} ${{PURGE_DUPS_BED}} ${{REFERENCE}};"
        " for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done"


rule merge_pri_hapdups_with_alt: #
    input:
        reference=out_dir_path  / ("contig/{assembler}/%s.contig.{assembler}.hap2.fasta" % config["genome_name"]),
        pri_hapdups=out_dir_path / ("{assembly_stage}/{assembler}/hap1/%s.{assembly_stage}.{assembler}.hap1.hap.fasta" % config["genome_name"])
    output:
        alt_plus_pri_hapdup=out_dir_path  / ("{assembly_stage}/{assembler}/input/%s.contig.{assembler}.hap2.fasta" % config["genome_name"]),
    log:
        std=output_dict["log"]  / "merge_pri_hapdups_with_alt.{assembler}.{assembly_stage}.log",
        cluster_log=output_dict["cluster_log"] / "merge_pri_hapdups_with_alt.{assembler}.{assembly_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_pri_hapdups_with_alt.{assembler}.{assembly_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt.{assembler}.{assembly_stage}..benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["merge_pri_hapdups_with_alt"] ,
        time=parameters["time"]["merge_pri_hapdups_with_alt"],
        mem=parameters["memory_mb"]["merge_pri_hapdups_with_alt"]
    threads: parameters["threads"]["merge_pri_hapdups_with_alt"]

    shell:
        " cat {input.reference} {input.pri_hapdups} > {output.alt_plus_pri_hapdup}"

rule create_link_for_purged_fasta:
    input:
        purged=rules.purge_dups.output.purged
    output:
        purged=out_dir_path  / ("{assembly_stage}/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"])
    log:
        std=output_dict["log"]  / "create_link_for_purged_fasta.{assembler}.{assembly_stage}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "create_link_for_purged_fasta.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_link_for_purged_fasta.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln {input.purged} {output.purged} > {log.std} 2>&1;"

