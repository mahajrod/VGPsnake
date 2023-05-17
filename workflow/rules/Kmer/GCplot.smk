rule gc_count:
    input:
        db=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"
    output:
        counts=output_dict["kmer"] / "{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.counts",
    #params:
        #kmer_length=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["kmer_length"],
        #min_coverage=lambda wildcards: "greater-than {0}".format(parameters["tool_options"]["gcp"][wildcards.datatype]["min_coverage"]) if "min_coverage" in parameters["tool_options"]["gcp"][wildcards.datatype] else "",
        #max_coverage=lambda wildcards: "less-than {0}".format(parameters["tool_options"]["gcp"][wildcards.datatype]["max_coverage"]) if "max_coverage" in parameters["tool_options"]["gcp"][wildcards.datatype] else "",
    log:
        gc_count=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.log",
        meryl=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.meryl.log",
        sort=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.sort.log",
        uniq=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.uniq.log",
        sed=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.sed.log",
        awk=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.awk.log",
        cluster_log=output_dict["cluster_log"] / "gc_plot.{datatype}.{stage}.L{min_coverage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gc_plot{datatype}.{stage}.L{min_coverage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["gc_count"],
        time=parameters["time"]["gc_count"],
        mem=parameters["memory_mb"]["gc_count"],
    threads:
        parameters["threads"]["gc_count"]
    shell: # output: coverage\tgc\tcount\n
         " meryl threads={threads} memory={resources.mem}m greater-than {wildcards.min_coverage} "
         " print {input.db} 2>{log.meryl} | count_kmer_gc.py 2>{log.gc_count} | sort -k2,2n -k1,1n 2>{log.sort} | "
         " uniq -c 2>{log.uniq} |  sed 's/^\s\+//;s/ /\t/' 2>{log.sed} | "
         " awk '{{printf \"%i\t%i\t%i\n\", $3,$2,$1 }}' > {output.counts} 2>{log.awk} "
"""
rule gc_plot:
    input:
        db=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"
    output:
        png=output_dict["kmer"] / "{datatype}/{stage}/kat/{datatype}.{stage}.{kmer_length}.jellyfish.kat.gcp.mx.png",
    params:
        #kmer_length=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["kmer_length"],
        hash_size=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["hash_size"],
        min_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["min_coverage"],
        max_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["max_coverage"],
        increment=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["increment"]
    log:
        gc_count=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.log",
        meryl=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}..meryllog",
        cluster_log=output_dict["cluster_log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gc_plot{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gc_plot.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kat"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kat"]["yaml"])
    resources:
        cpus=parameters["threads"]["gc_plot"],
        time=parameters["time"]["gc_plot"],
        mem=parameters["memory_mb"]["gc_plot"],
    threads:
        parameters["threads"]["kat_gcp"]
    shell:
         " KMER_LEN={wildcards.kmer_length}; "
         " meryl threads={threads} memory={resources.mem}m print {input.db} 2>{log.meryl} | "
         " count_kmer_gc.py | sort -k1,1n -k2,2n |  uniq -c |  sed 's/^\s\+//;s/ /\t/' | "
         " awk '{{printf \"%i\t%i\t%i\n\", $2,$3,$1 }}' > {output.kmer} 2>{log.gc_count};"
"""