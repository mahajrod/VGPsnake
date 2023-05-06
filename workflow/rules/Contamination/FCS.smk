
rule fcs: #
    priority: 1000
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        db=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["path"],
        image=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["image_path"],
    output:
        report=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.taxonomy",
        summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.summary"
        #report=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.taxonomy.txt",
        #summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.fcs_gx_report.txt"
    params:
        tax_id=config["tax_id"]
    log:
        std=output_dict["log"]  / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.log",
        cluster_log=output_dict["cluster_log"] / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.err"
    benchmark:
        output_dict["benchmark"]  / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        cpus=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["threads"],
        time=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["time"],
        mem=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["memory_mb"],
        fcs=lambda wildcards: 1 if wildcards.database == "fcs_gx" else 0
    threads: lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["threads"],

    shell:
        " OUTDIR=`dirname {output.report}`; "
        " export SINGULARITYENV_TMPDIR=${{OUTDIR}}; "
        " export SINGULARITYENV_SQLITE_TMPDIR=${{OUTDIR}}; "
        " NUM_CORES={threads}; "
        " export FCS_DEFAULT_IMAGE={input.image}; "
        " fcs.py  screen genome --fasta {input.fasta} --out-dir `dirname {output.report}` --tax-id {params.tax_id} --gx-db {input.db} > {log.std} 2>&1; "
        " REPORT={output.report}; "
        " SUMMARY={output.summary}; "
        " mv ${{REPORT%.{wildcards.database}.taxonomy}}.{params.tax_id}.{wildcards.database}_report.txt ${{REPORT}}; "
        " mv ${{SUMMARY%.{wildcards.database}.summary}}.{params.tax_id}.taxonomy.rpt ${{SUMMARY}}; "

rule fcs_adaptor: #
    priority: 2000
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        #db=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["path"],
        image=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["image_path"],
    output:
        report=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.report",
        report_jsonl=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.report.jsonl",
        #summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.summary"
        #report=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.taxonomy.txt",
        #summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.fcs_gx_report.txt"
    params:
        tax_id=config["tax_id"]
    log:
        std=output_dict["log"]  / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.log",
        cluster_log=output_dict["cluster_log"] / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.err"
    benchmark:
        output_dict["benchmark"]  / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        cpus=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["threads"],
        time=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["time"],
        mem=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["memory_mb"],
    threads: lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["threads"],

    shell:
        " OUTDIR=`dirname {output.report}`; "
        " export SINGULARITYENV_TMPDIR=${{OUTDIR}}; "
        " export SINGULARITYENV_SQLITE_TMPDIR=${{OUTDIR}}; "
        " run_fcsadaptor.sh --image {input.image} --fasta-input {input.fasta} --output-dir `dirname {output.report}` --prok --container-engine singularity > {log.std} 2>&1; "
        " mv `dirname {output.report}`/fcs_adaptor_report.txt {output.report}; "
        " mv `dirname {output.report_jsonl}`/combined.calls.jsonl {output.report_jsonl}; "
