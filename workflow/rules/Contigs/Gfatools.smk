rule gfa2fasta:
    input:
        gfa=output_dict["contig"] / ("{assembler}/%s.contig.gfa2fasta.pacbio.hic.{haplotype}_ctg.gfa" % config["genome_name"])
    output:
        fasta=output_dict["contig"] / ("{assembler}/%s.contig.gfa2fasta.pacbio.hic.{haplotype}_ctg.fasta" % config["genome_name"])

    log:
        std=output_dict["log"] / "gfa2fasta.{assembler}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "gfa2fasta.cluster.{assembler}.{haplotype}.log",
        cluster_err=output_dict["cluster_error"] / "gfa2fasta.cluster.{assembler}.{haplotype}.err"
    benchmark:
        output_dict["benchmark"] / "gfa2fasta.{assembler}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["gfa2fasta"],
        time=parameters["time"]["gfa2fasta"],
        mem=parameters["memory_mb"]["gfa2fasta"],
    threads:
        parameters["threads"]["gfa2fasta"]
    shell:
         " gfatools gfa2fa {input.gfa} > {output.fasta} 2>{log.std};"

