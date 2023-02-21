rule meryl:
    input:
        output_dict["data"] / ("fastq/{datatype}}/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{fileprefix}.{kmer_length}.meryl"),
        #histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.histo"
    params:
        kmer_length=lambda wildcards: parameters["tool_options"]["meryl"][wildcards.datatype]["kmer_length"],
        hash_size=lambda wildcards: parameters["tool_options"]["meryl"][wildcards.datatype]["hash_size"],
        min_coverage=lambda wildcards: parameters["tool_options"]["meryl"][wildcards.datatype]["min_coverage"],
        max_coverage=lambda wildcards: parameters["tool_options"]["meryl"][wildcards.datatype]["max_coverage"],
        increment=lambda wildcards: parameters["tool_options"]["meryl"][wildcards.datatype]["increment"]
    log:
        std=output_dict["log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=parameters["memory_mb"]["meryl"],
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl k={params.kmer_length} cpu={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"


rule merge_meryl:
    input:
        lambda wildcards:
            expand(output_dict["kmer"] / ("%s/%s/%s.%s.{fileprefix}.%i.meryl" % (wildcards.datatype,
                                                                                 wildcards.stage,
                                                                                 wildcards.datatype,
                                                                                 wildcards.stage,
                                                                                 parameters["tool_options"]["meryl"][wildcards.datatype]["kmer_length"],)),
                   fileprefix=input_file_prefix_dict[wildcards.datatype],
                   allow_missing=True,)
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"),
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.histo"

    log:
        count_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.count.log",
        histo_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.histo.log",
        cluster_log=output_dict["cluster_log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=parameters["memory_mb"]["meryl"],
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl cpu={threads} memory={resources.mem}m count "
         " union-sum output {output.db_dir} {input} 1>{log.count_log} 2>&1;"
         " meryl cpu={threads} memory={resources.mem}m "
         " histogram output {output.histo} {output.db_dir} 1>{log.histo_log} 2>&1"

