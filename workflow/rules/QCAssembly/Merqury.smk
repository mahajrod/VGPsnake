rule merqury:
    input:
        meryl_db_dir=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{fileprefix}.{kmer_length}.meryl",
        primary_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.p_ctg.fasta" % config["genome_name"]),
        alternative_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.a_ctg.fasta" % config["genome_name"])
    output:
        #summary=output_dict["assembly_qc"] /("{assembly_stage}/{assembler}/{assembly_stage}.{assembler}.pacbio.hic.{haplotype}_ctg.gfa" % config["genome_name"]),
        dir=directory(output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}"
                                                   % config["genome_name"]))
    params:
        out_prefix=
    log:
        std=output_dict["log"] / "merqury.{assembler}.{assembly_stage}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "merqury.{assembler}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merqury.{assembler}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merqury.{assembler}.{assembly_stage}.{haplotype}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["merqury"],
        time=parameters["time"]["merqury"],
        mem=parameters["memory_mb"]["merqury"],
    threads:
        parameters["threads"]["merqury"]
    shell:
         " OMP_NUM_THREADS={threads} merqury.sh {input.meryl_db_dir} "
         " {input.primary_assembly} {input.alternative_assembly} {params.out_prefix}  1>{log.std} 2>&1;"
         
         
ln -s $MERQURY/merqury.sh		# Link merqury
./merqury.sh <read-db.meryl> [<mat.meryl> <pat.meryl>] <asm1.fasta> [asm2.fasta] <out>

Usage: merqury.sh <read-db.meryl> [<mat.meryl> <pat.meryl>] <asm1.fasta> [asm2.fasta] <out>
	<read-db.meryl>	: k-mer counts of the read set
	<mat.meryl>		: k-mer counts of the maternal haplotype (ex. mat.only.meryl or mat.hapmer.meryl)
	<pat.meryl>		: k-mer counts of the paternal haplotype (ex. pat.only.meryl or pat.hapmer.meryl)
	<asm1.fasta>	: Assembly fasta file (ex. pri.fasta, hap1.fasta or maternal.fasta)
	[asm2.fasta]	: Additional fasta file (ex. alt.fasta, hap2.fasta or paternal.fasta)
	*asm1.meryl and asm2.meryl will be generated. Avoid using the same names as the hap-mer dbs
	<out>		: Output prefix
