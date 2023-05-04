rule meryl_assembly:
    input:
        fasta=out_dir_path / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}.fasta"
    output:
        db_dir=directory(out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}")
    log:
        std=output_dict["log"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl_assembly"],
        time=parameters["time"]["meryl_assembly"],
        mem=parameters["memory_mb"]["meryl_assembly"],
    threads:
        parameters["threads"]["meryl_assembly"]
    shell:
         " meryl k={wildcards.assembly_kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"

#def get_rest_haplotype_list()

rule meryl_extract_unique_hap_kmers:
    input:
        target_hap_db_dir=out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}",
        #hap1_db_dir=out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.hap1.{assembly_kmer_length}",
        rest_hap_db_dirs=lambda wildcards: expand(out_dir_path / "{0}/{1}/kmer/{2}.{0}.{3}.{4}".format(wildcards.stage,
                                                                                                       wildcards.parameters,
                                                                                                       wildcards.genome_prefix,
                                                                                                       wildcards.haplotype,
                                                                                                       wildcards.assembly_kmer_length) ,
                                                 haplotype=set(haplotype_list) - set(wildcards.haplotype),
                                                 allow_missing=True)
    output:
        unique_hap_db_dir=out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.unique",
        #unique_hap2_db_dir=out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.hap2.{assembly_kmer_length}.unique"
    log:
        std=output_dict["log"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.log",
        #hap2=output_dict["log"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.hap2.log",
        cluster_log=output_dict["cluster_log"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl_extract_unique_hap_kmers"],
        time=parameters["time"]["meryl_extract_unique_hap_kmers"],
        mem=parameters["memory_mb"]["meryl_extract_unique_hap_kmers"],
    threads:
        parameters["threads"]["meryl_extract_unique_hap_kmers"]
    shell:
         " meryl threads={threads} memory={resources.mem}m difference {input.target_hap_db_dir} {input.rest_hap_db_dirs} output {output.unique_hap_db_dir} > {log.std} 2>&1; "
         #" meryl threads={threads} memory={resources.mem}m difference {input.hap2_db_dir} {input.hap1_db_dir} output {output.unique_hap2_db_dir} > {log.hap2} 2>&1; "


rule extract_reads_by_unique_hap_kmers:
    input:
        rest_hap_db_dirs=lambda wildcards: expand(out_dir_path / "{0}/{1}/kmer/{2}.{0}.{3}.{4}.unique".format(wildcards.stage,
                                                                                                              wildcards.parameters,
                                                                                                              wildcards.genome_prefix,
                                                                                                              wildcards.haplotype,
                                                                                                              wildcards.assembly_kmer_length) ,
                                                 haplotype=set(haplotype_list) - set(wildcards.haplotype),
                                                 allow_missing=True),
        forward_read=output_dict["data"]  / ("fastq/hic/raw/{pairprefix}%s%s" % (input_forward_suffix_dict["hic"], config["fastq_extension"])),
        reverse_read=output_dict["data"]  / ("fastq/hic/raw/{pairprefix}%s%s" % (input_reverse_suffix_dict["hic"], config["fastq_extension"])),
    output:
        forward_hap1_read=out_dir_path / "{stage}/{parameters}/fastq/{haplotype}/{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.{haplotype}_1.fastq.gz",
        reverse_hap1_read=out_dir_path / "{stage}/{parameters}/fastq/{haplotype}/{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.{haplotype}_2.fastq.gz",
        #forward_hap2_read=out_dir_path / "{stage}/{parameters}/kmer/{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.hap2_1.fastq.gz",
        #reverse_hap2_read=out_dir_path / "{stage}/{parameters}/kmer/{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.hap2_2.fastq.gz",
    log:
        std=output_dict["log"] / "extract_reads_by_unique_hap_kmers.{stage}.{parameters}.{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.{haplotype}.log",
        #hap2=output_dict["log"] / "extract_reads_by_unique_hap_kmers.{stage}.{parameters}.{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.hap2.log",
        cluster_log=output_dict["cluster_log"] / "extract_reads_by_unique_hap_kmers.{stage}.{parameters}.{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_reads_by_unique_hap_kmers.{stage}.{parameters}.{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "extract_reads_by_unique_hap_kmers.{stage}.{parameters}.{pairprefix}.{genome_prefix}.AK{assembly_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_reads_by_unique_hap_kmers"],
        time=parameters["time"]["extract_reads_by_unique_hap_kmers"],
        mem=parameters["memory_mb"]["extract_reads_by_unique_hap_kmers"],
    threads:
        parameters["threads"]["extract_reads_by_unique_hap_kmers"]
    shell:
         " meryl-lookup -exclude -sequence {input.forward_read} {input.reverse_read} "
         " -mers {input.rest_hap_db_dirs} -output {output.forward_hap1_read} {output.reverse_hap1_read} > {log.std} 2>&1;"
         # " meryl-lookup -exclude -sequence {input.forward_read} {input.reverse_read} "
         #" -mers {input.unique_hap2_db_dir} -output {output.forward_hap1_read} {output.reverse_hap1_read} > {log.hap1} 2>&1;"

