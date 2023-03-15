
rule pretextmap: #
    input:
        bam=out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/alignment/%s.hic_scaffolding.{assembler}.{haplotype}.bwa.filtered.rmdup.bam"  % config["genome_name"]),
    output:
        map=out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/alignment/%s.hic_scaffolding.{assembler}.{haplotype}.bwa.filtered.rmdup.map.pretext"  % config["genome_name"]),
    params:
        min_mapq=parameters["tool_options"]["pretextmap"]["mapq"],
        sortby=parameters["tool_options"]["pretextmap"]["sortby"],
        sortorder=parameters["tool_options"]["pretextmap"]["sortorder"]
    log:
        view=output_dict["log"]  / "pretextmap.{assembler}.hic_scaffolding.{haplotype}.view.log",
        map=output_dict["log"]  / "pretextmap.{assembler}.hic_scaffolding.{haplotype}.map.log",
        cluster_log=output_dict["cluster_log"] / "pretextmap.{assembler}.hic_scaffolding.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextmap.{assembler}.hic_scaffolding.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextmap.{assembler}.hic_scaffolding.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["pretextmap"] ,
        time=parameters["time"]["pretextmap"],
        mem=parameters["memory_mb"]["pretextmap"]
    threads: parameters["threads"]["pretextmap"]

    shell:
        " samtools view -h {input.bam} 2>{log.view} | "
        " PretextMap -o {output.map} --sortby {params.sortby} --sortorder {params.sortorder} "
        " --mapq {params.min_mapq} > {log.map} 2>&1"

rule pretextsnapshot: 
    input:
        map=rules.pretextmap.output.map
    output:
        dir=directory(out_dir_path / ("hic_scaffolding/{assembler}/{haplotype}/alignment/%s.map.{assembler}.{haplotype}.{resolution}.{ext}" % config["genome_name"])),
        #fig=out_dir_path / ("hic_scaffolding/{assembler}/{haplotype}/alignment/%s.hic_scaffolding.{assembler}.{haplotype}.bwa.filtered.rmdup.map.{resolution}.{ext}" % config["genome_name"]),
    params:
        #dir=lambda wildcards: out_dir_path / "hic_scaffolding/{0}/{1}/alignment/{2}.map.{3}.{4}".format(config["genome_name"],
        #                                                                                                wildcards.assembler,
        #                                                                                                wildcards.haplotype,
        #                                                                                                wildcards.resolution,
        #                                                                                                wildcards.ext),
        #prefix=lambda wildcards: "{0}.hic_scaffolding.{1}.{2}.bwa.filtered.rmdup.map.{3}".format(config["genome_name"],
        #                                                                                         wildcards.assembler,
        #                                                                                         wildcards.haplotype,
        #                                                                                         wildcards.resolution),
        sequences=parameters["tool_options"]["pretextsnapshot"]["sequences"],
    log:
        std=output_dict["log"]  / "pretextsnapshot.{assembler}.hic_scaffolding.{haplotype}.{ext}.{resolution}.log",
        cluster_log=output_dict["cluster_log"] / "pretextsnapshot.{assembler}.hic_scaffolding.{haplotype}.{ext}.{resolution}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextsnapshot{assembler}.hic_scaffolding.{haplotype}.{ext}.{resolution}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextsnapshot.{assembler}.hic_scaffolding.{haplotype}.{ext}.{resolution}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["pretextmap"] ,
        time=parameters["time"]["pretextmap"],
        mem=parameters["memory_mb"]["pretextmap"]
    threads: parameters["threads"]["pretextmap"]
    shell:
        " PretextSnapshot --sequences {params.sequences} -r {wildcards.resolution} -f {wildcards.ext} "
        " -m {input.map} -o {output.dir}  > {log.std} 2>&1"

