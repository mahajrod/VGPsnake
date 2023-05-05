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
        cpus=parameters["threads"]["smudgeplot_plot"],
        time=parameters["time"]["smudgeplot_plot"],
        mem=parameters["memory_mb"]["smudgeplot_plot"],
    threads:
        parameters["threads"]["smudgeplot_plot"]
    shell:
         " LOWER_BOUNDARY=`smudgeplot.py cutoff {input.histo} L 2>{log.lower}`; "
         " UPPER_BOUNDARY=`smudgeplot.py cutoff {input.histo} U 2>{log.upper}`; "
         " echo -e \"low_boundary\tupper_boundary\n${{LOWER_BOUNDARY}}\t${{UPPER_BOUNDARY}}\n\" > {output.boundaries}"

rule smudgeplot_hetkmers:
    input:
        kmer=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer",
    output:
        coverages=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_coverages.tsv",
        sequences=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_sequences.tsv",
        #smudgeplot=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_smudgeplot.png",
        #summary=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_summary_table.tsv"
    log:
        hetkmers=output_dict["log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.hetkmers.log",
        cluster_log=output_dict["cluster_log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["smudgeplot_hetkmers"],
        time=parameters["time"]["smudgeplot_hetkmers"],
        mem=parameters["memory_mb"]["smudgeplot_hetkmers"],
        smudgeplot_hetkmers=1
    threads:
        parameters["threads"]["smudgeplot_hetkmers"]
    shell:
         " COV_OUT={output.coverages}; "
         " PREFIX=${{COV_OUT%_coverages.tsv}}; "
         #" HAPLOID_COVERAGE=`awk 'NR==2 {{print 2 * $2}}' {input.genomescope_report}`; "
         " smudgeplot.py hetkmers -o ${{PREFIX}} {input.kmer} > {log.hetkmers} 2>&1; "
         #" smudgeplot.py plot -k {wildcards.kmer_length} -n ${{HAPLOID_COVERAGE}}  -o ${{PREFIX}} {output.coverages} > {log.plot} 2>&1; "

rule smudgeplot_plot:
    input:
        coverages=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_coverages.tsv",
        genomescope_report=output_dict["kmer"] / ("{datatype}/{stage}/genomescope/%s.%s.filtered.%s.%s.genomescope.parameters" % (config["genome_prefix"],
                                                                                                                                  config["final_kmer_datatype"],
                                                                                                                                  config["final_kmer_length"],
                                                                                                                                  config["final_kmer_counter"])),
    output:
        smudgeplot=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_smudgeplot.png",
        summary=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_summary_table.tsv",
        smudgeplot_no_priors=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.no_priors_smudgeplot.png",
        summary_no_priors=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.no_priors_summary_table.tsv"
    log:
        plot=output_dict["log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.plot.log",
        cluster_log=output_dict["cluster_log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["smudgeplot_plot"],
        time=parameters["time"]["smudgeplot_plot"],
        mem=parameters["memory_mb"]["smudgeplot_plot"],
    threads:
        parameters["threads"]["smudgeplot_plot"]
    shell:
         " COV_OUT={output.smudgeplot}; "
         " PREFIX=${{COV_OUT%_smudgeplot.png}}; "
         " HAPLOID_COVERAGE=`awk 'NR==2 {{print $2}}' {input.genomescope_report}`; "
         #" smudgeplot.py hetkmers -o ${{PREFIX}} {input.kmer} > {log.hetkmers} 2>&1; "
         " smudgeplot.py plot -k {wildcards.kmer_length} -n ${{HAPLOID_COVERAGE}}  -o ${{PREFIX}} {input.coverages} > {log.plot} 2>&1; "
         " COV_OUT_NO_PRIORS={output.smudgeplot_no_priors}; "
         " PREFIX_NO_PRIORS=${{COV_OUT_NO_PRIORS%_smudgeplot.png}}; "
         " smudgeplot.py plot -k {wildcards.kmer_length} -o ${{PREFIX_NO_PRIORS}} {input.coverages} > {log.plot} 2>&1; "


rule compress_kmer:
    input:
        kmer=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer",
        summary=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_summary_table.tsv"
    output:
        kmer_gz=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer.gz"
    log:
        std=output_dict["log"] / "compress_kmer.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.log",
        cluster_log=output_dict["cluster_log"] / "compress_kmer{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "compress_kmer.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "compress_kmer.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["compress_kmer"],
        time=parameters["time"]["compress_kmer"],
        mem=parameters["memory_mb"]["compress_kmer"],
        smudgeplot=1
    threads:
        parameters["threads"]["compress_kmer"]
    shell:
         " pigz -p 10 {input.kmer} > {log.std} 2>&1; "
