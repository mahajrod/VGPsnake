localrules: parse_genomescope_output

rule genomescope:
    input:
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.histo"
    output:
        dir=directory(output_dict["kmer"] / "{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}"),
    params:
        ploidy=config["ploidy"],
        genome_name=config["genome_name"],
        max_coverage=lambda wildcards: parameters["tool_options"][wildcards.kmer_tool][wildcards.datatype]["max_coverage"],

    log:
        std=output_dict["log"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.count.log",
        cluster_log=output_dict["cluster_log"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["genomescope"],
        time=parameters["time"]["genomescope"],
        mem=parameters["memory_mb"]["genomescope"],
    threads:
        parameters["threads"]["genomescope"]
    shell:
         " genomescope.R  -i {input.histo} -p {params.ploidy} -k {wildcards.kmer_length}  "
         " -n {params.genome_name}  -m {params.max_coverage}  --fitted_hist  --testing  -o {output.dir} > log.std 2>&1"


rule parse_genomescope_output:
    input:
        summary=output_dict["kmer"] / ("{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/%s_summary.txt" % config["genome_name"]),
        model=output_dict["kmer"] / ("{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/%s_model.txt" % config["genome_name"])
    output:
        output_dict["kmer"] / "{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters"
    log:
        std=output_dict["log"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.count.log",
        cluster_log=output_dict["cluster_log"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["parse_genomescope_output"],
        time=parameters["time"]["parse_genomescope_output"],
        mem=parameters["memory_mb"]["parse_genomescope_output"],
    threads:
        parameters["threads"]["parse_genomescope_output"]
    shell:
         " GENLEN=`grep 'Genome Haploid Length' {input.summary} | sed 's/,//g;s/ \\{2,\\}/\t/g' | cut -f 3 | sed 's/ .*//'`   "
         " LAMBDA=`grep 'kmercov' {input.model} | tail -n  1 | awk '{printf \"%.0f\", \$2}`"
         " echo -e \"Genome size\t${GENLEN}\nLambda\t${LAMBDA}\n\" > {output}"
