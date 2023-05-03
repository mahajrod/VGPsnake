localrules: smudgeplot_assess

rule smudgeplot_assess:
    input:
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.histo"
    output:
        boundaries=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.smudgeplot.boundaries",
    log:
        upper=output_dict["log"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.upper.log",
        lower=output_dict["log"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.lower.log",
        cluster_log=output_dict["cluster_log"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["smudgeplot"],
        time=parameters["time"]["smudgeplot"],
        mem=parameters["memory_mb"]["smudgeplot"],
    threads:
        parameters["threads"]["smudgeplot"]
    shell:
         " LOWER_BOUNDARY=`smudgeplot.py cutoff {input.histo} L 2>{log.lower}`; "
         " UPPER_BOUNDARY=`smudgeplot.py cutoff {input.histo} U 2>{log.upper}`; "
         " echo -e \"low_boundary\tupper_boundary\n${{LOWER_BOUNDARY}}\t${{UPPER_BOUNDARY}}\n\" > {output.boundaries}"

rule smudgeplot:
    input:
        kmer=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.kmer.gz",
        genomescope_report=output_dict["kmer"] / ("{datatype}/{stage}/genomescope/%s.%s.filtered.%s.%s.genomescope.parameters" % (config["genome_prefix"],
                                                                                                                                  config["final_kmer_datatype"],
                                                                                                                                  config["final_kmer_length"],
                                                                                                                                  config["final_kmer_counter"])),
    output:
        coverages=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_coverages.tsv",
        sequences=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_sequences.tsv",
        smudgeplot=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_smudgeplot.png",
        summary=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_summary_table.tsv"
    log:
        hetkmers=output_dict["log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.hetkmers.log",
        plot=output_dict["log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.plot.log",
        cluster_log=output_dict["cluster_log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["smudgeplot"],
        time=parameters["time"]["smudgeplot"],
        mem=parameters["memory_mb"]["smudgeplot"],
    threads:
        parameters["threads"]["smudgeplot"]
    shell:
         " COV_OUT={output.coverages}; "
         " PREFIX=${{COV_OUT%_coverages.tsv}}; "
         " HAPLOID_COVERAGE=`awk 'NR==2 {{print 2 * $2}}'`; "
         " smudgeplot.py hetkmers -o ${{PREFIX}} <(zcat {input.kmer}) > {log.hetkmers} 2>&1; "
         " smudgeplot.py plot -k {wildcards.kmer_length} -n ${{HAPLOID_COVERAGE}}  -o ${{PREFIX}} {output.coverages} > {log.plot} 2>&1; "
