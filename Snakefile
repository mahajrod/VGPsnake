import os
#import logging
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath

import yaml

import pandas as pd

#logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

#---- Functions ----
def p_distance(seq_a, seq_b, seq_len):
    dist = 0
    for i in range(0, seq_len):
        if seq_a[i] != seq_b[i]:
            dist += 1
    return dist

def get_common_prefix_ans_suffixes(seq_a, seq_b):
    seq_a_len = len(seq_a)
    seq_b_len = len(seq_b)

    prefix = ""
    for i in range(0, min(seq_a_len, seq_b_len)):
        if seq_a[i] != seq_b[i]:
           return prefix, seq_a[i:], seq_b[i:]
        prefix += seq_a[i]
    return prefix, "", ""

def convert_posixpath2str_in_dict(dictionary):
    output_dictionary = deepcopy(dictionary)
    for entry in output_dictionary:
        if isinstance(output_dictionary[entry], PosixPath):
            output_dictionary[entry] = str(output_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            output_dictionary[entry] = convert_posixpath2str_in_dict(output_dictionary[entry])

    return output_dictionary

def find_cmap(bionano_dir, cmap_extension=".cmap"): # TODO: modify when input for bionano will be clear
    bionano_dir_path = bionano_dir if isinstance(bionano_dir, PosixPath) else Path(bionano_dir)
    cmap_list = list(bionano_dir_path.glob("*{0}".format(cmap_extension)))
    if len(cmap_list) > 1:
        raise ValueError(
                         "ERROR!!! More than one cmap file was found: {0}".format(
                                                                                  ", ".join(list(map(str, cmap_list)))
                                                                                  )
                         )
    return cmap_list[0]

def find_fastqs(fastq_dir, fastq_extension=".fastq.gz"):
    fastq_dir_path = fastq_dir if isinstance(fastq_dir, PosixPath) else Path(fastq_dir)
    return  sorted(list(fastq_dir_path.glob("*{0}".format(fastq_extension))))


def copy_absent_entries(input_dictionary, output_dictionary):
    for entry in input_dictionary:
        if entry not in output_dictionary:
            output_dictionary[entry] = deepcopy(input_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            copy_absent_entries(input_dictionary[entry], output_dictionary[entry])

#----


#-- Initialization of path variables from config file --
#logging.info("Initialization of path variables...")
#---- Initialization of path variables for input----
input_dir_path = Path(config["input_dir"])

input_dict = {}
data_types = config["data_types"].split(",")

for datatype in data_types:
    input_dict[datatype] = {}
    input_dict[datatype]["dir"] = input_dir_path / datatype
    input_dict[datatype]["run_dir"] = input_dict[datatype]["dir"] / "run"
    if datatype == "bionano":
        input_dict[datatype]["cmap"] = None # TODO: implement parsing bionano .cmap filename
    else:
        input_dict[datatype]["fastq_dir"] = input_dict[datatype]["dir"] / "fastq"

#----
#---- Initialization of path variables for output ----
out_dir_path = Path(config["out_dir"])
output_dict = {}

for first_level_sub_dir in config["first_level_subdir_list"]:
    output_dict[first_level_sub_dir] = out_dir_path / config["{0}_subdir".format(first_level_sub_dir)]

#----
#---- Initialization path variables for resources ----
#----
#---- Setting mode of pipeline ----
#.info("Setting and adjusting pipeline mode...")

#pipeline_mode = config["mode"]
#starting_point = config["starting_point"]

#-------- Verification of input datatypes --------

fastq_based_data_type_set = set(data_types) & set(config["fastq_based_data"])
genome_size_estimation_data_type_set = set(config["genome_size_estimation_data"]) & fastq_based_data_type_set

#logging.info("Verifying datatypes...")
for d_type in data_types:
    if d_type not in config["allowed_data_types"]:
        #logging.error("Unknown data type: {0}".format(d_type))
        raise ValueError("ERROR!!! Unknown data type: {0}".format(d_type))

#--------

#----

#---- Checking input files ----
#logging.info("Checking input files...")

input_filedict = {}
input_file_prefix_dict = {}
input_forward_suffix_dict = {}
input_reverse_suffix_dict = {}
input_pairprefix_dict = {}

for d_type in fastq_based_data_type_set:
    input_filedict[d_type] = find_fastqs(input_dict[d_type]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict[d_type] = list(map(lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                                input_filedict[d_type]))
# check filenames of paired data
for d_type in set(config["paired_fastq_based_data"]) & fastq_based_data_type_set:
   if (len(input_filedict[d_type]) % 2) != 0:
        raise ValueError("ERROR!!! {0} fastq files seems to be unpaired or misrecognized".format(d_type))
   for forward, reverse in zip(input_filedict[d_type][::2], input_filedict[d_type][1::2]):
        #print(forward, reverse)
        if p_distance(str(forward), str(reverse), len(str(forward))) > 1:
            raise ValueError("ERROR!!! Forward and reverse read files differs by more than one symbol:\n\t{0}\n\t{1}".format(str(forward),
                                                                                                                             str(reverse)))
#get_suffixes for paired fastq data
for d_type in set(config["paired_fastq_based_data"]) & fastq_based_data_type_set:
    input_forward_suffix_dict[d_type] = set()
    input_reverse_suffix_dict[d_type] = set()
    input_pairprefix_dict[d_type] = []
    for forward_prefix, reverse_prefix in zip(input_file_prefix_dict[d_type][::2], input_file_prefix_dict[d_type][1::2]):
        common_prefix, forward_suffix, reverse_suffix = get_common_prefix_ans_suffixes(forward_prefix, reverse_prefix)
        input_pairprefix_dict[d_type].append(common_prefix)
        input_forward_suffix_dict[d_type].add(forward_suffix)
        input_reverse_suffix_dict[d_type].add(reverse_suffix)
    if (len(input_forward_suffix_dict[d_type]) > 1) or (len(input_reverse_suffix_dict[d_type]) > 1):
        raise ValueError("ERROR!!! Multiple different suffixes in filenames of %s data!" % d_type)
    input_forward_suffix_dict[d_type] = list(input_forward_suffix_dict[d_type])[0]
    input_reverse_suffix_dict[d_type] = list(input_reverse_suffix_dict[d_type])[0]


if "bionano" in data_types: # TODO: modify when input for bionano will be clear
    input_filedict["bionano"] = find_cmap(input_dict["bionano"]["dir"], cmap_extension=config["cmap_extension"])


#---- Initialize tool parameters ----
#logging.info("Initializing tool parameters...")
#check if custom restriction sites were provided:
if config["custom_enzyme_set"] is not None:
    config["parameters"]["default"]["tool_options"]["salsa2"]["restriction_seq"]["custom"] = config["custom_enzyme_set"]
    if "tool_options" in config["parameters"][config["parameter_set"]]:
        if "salsa2" in config["parameters"][config["parameter_set"]]["tool_options"]:
            if "restriction_seq" in config["parameters"][config["parameter_set"]]["tool_options"]["salsa2"]:
                if "custom" in config["parameters"][config["parameter_set"]]["tool_options"]["salsa2"]["restriction_seq"]:
                    config["parameters"][config["parameter_set"]]["tool_options"]["salsa2"]["restriction_seq"]["custom"] = config["custom_enzyme_set"]
    config["hic_enzyme_set"] = "custom"

if config["parameter_set"] not in config["parameters"]:
    raise ValueError("Error!!! Unknown set of tool parameters: {0}".format(config["parameter_set"]))

copy_absent_entries(config["parameters"]["default"], config["parameters"][config["parameter_set"]]) # set default values for options absent in  "parameter_set"

for key in list(config["parameters"].keys()): # remove unused sets of parameters
    if key != config["parameter_set"]:
        config["parameters"].pop(key)

parameters = config["parameters"][config["parameter_set"]] # short alias for used set of parameters

#check if final_kmer_tool is present in "kmer_counter_list"
if config["final_kmer_counter"] not in config["kmer_counter_list"]:
    config["kmer_counter_list"].append(config["final_kmer_counter"])
    #logging.info("Warning! final_kmer_counter is not in kmer_counter_list! Added...")

#check if final_kmer_length is present in parameters of final_kmer_tool
for dat_type in genome_size_estimation_data_type_set:
    if config["final_kmer_length"] not in parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"]:
        parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"].append(config["final_kmer_length"])
        #logging.info("Warning! Final_kmer_length is not in parameters of final_kmer_counter! Added...")

#----
#---- Configure stages ----
config["stage_list"] = []

# Select configuration and combine stages from all mega_stages in a single list without nesting
if config["mode"] == "preprocessing":
    mega_stage_list = ["preprocessing"]
elif config["mode"] == "qc":
    mega_stage_list = ["preprocessing", "qc"]
elif config["mode"] == "assembly":
    mega_stage_list = ["preprocessing", "qc", "assembly"]
else:
    raise ValueError("ERROR!!! Unknown mode: %s" % config["mode"])

for mega_stage in mega_stage_list:
    custom_megastage_entry = "custom_" + mega_stage + "_stages"
    if (custom_megastage_entry in config) and (config[custom_megastage_entry]):
        config["stage_list"].append(config[custom_megastage_entry])
    else:
        config["stage_list"] += config["allowed_stage_list"][mega_stage][config[mega_stage + "_mode"]][config["starting_point"]]

stage_dict = OrderedDict()
for stage, stage_index in zip(config["stage_list"], range(0, len(config["stage_list"]))):
    stage_dict[stage] = OrderedDict()
    stage_dict[stage]["prev_stage"] = None if stage_index == 0 else config["stage_list"][stage_index-1]

#----


#---- Save configuration and input files ----
final_config_yaml = output_dict["config"] / "config.final.yaml"
final_input_yaml = output_dict["config"] / "input.final.yaml"

os.makedirs(output_dict["config"], exist_ok=True)

with open(final_config_yaml, 'w') as final_config_fd, open(final_input_yaml, 'w') as final_input_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False, sort_keys=False)
    yaml.dump(convert_posixpath2str_in_dict(input_dict), final_input_fd, default_flow_style=False, sort_keys=False)


#-------------------------------------------
localrules: all
ruleorder: create_fastq_links > fastqc

results_dict = {}

assembler_list = ["hifiasm", ] # TODO: implement possibility of other assemblers

haplotype_list = ["hap{0}".format(i) for i in range(1, config["ploidy"] + 1)]
primary_haplotype = "hap1"

results_list = []

#---- Create output filelist ----
if "check_reads" in config["stage_list"]:
    results_list += [
                     final_config_yaml,
                     final_input_yaml
                     ]

if "check_draft" in config["stage_list"]:
    results_list += [ ] # TODO: implement

if "read_qc" in config["stage_list"]:
    results_list += [*[expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=input_file_prefix_dict[dat_type], #list(
                                 #    map(
                                 #        lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                 #        input_filedict[dat_type])
                                 #        )
                               ) for dat_type in fastq_based_data_type_set],
                      expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                             datatype=fastq_based_data_type_set,
                             stage=["raw",]), ]

if "draft_qc" in config["stage_list"]:
    results_list += [ ] # TODO: implement

if "filter_reads" in config["stage_list"]:
    results_list += [expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                                    fileprefix=input_file_prefix_dict["hifi"]) if "hifi" in fastq_based_data_type_set else [],
                             expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                                    datatype=["hifi", ],
                                    stage=["filtered", ],
                                    fileprefix=input_file_prefix_dict["hifi"],
                                    ) if "hifi" in fastq_based_data_type_set else [],
                             expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                                    datatype=["hifi"],
                                    stage=["filtered",]) if "hifi" in fastq_based_data_type_set else [], # only hifi filtration was implemented yet
                             *[[expand(output_dict["kmer"] / "{datatype}/{stage}/genomescope/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters",
                                    datatype=[dat_type,],
                                    genome_prefix=[config["genome_prefix"], ],
                                    stage=["filtered",],
                                    kmer_tool=[kmer_tool,],
                                    kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                                    ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set]
                              ]

if "filter_draft" in config["stage_list"]:
    results_list += [ ] # TODO: implement

if "contig" in config["stage_list"]:
    assembler_list = config["stage_coretools"]["contig"][config["contig_datatype"]]
    stage_dict["contig"]["parameters"] = {}

    for assembler in assembler_list:
        for option_set in config["coretool_option_sets"][assembler]:
            parameters_label="{0}_{1}".format(assembler, option_set)
            stage_dict["contig"]["parameters"][parameters_label] = {}
            stage_dict["contig"]["parameters"][parameters_label]["assembler"] = assembler
            stage_dict["contig"]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][assembler][option_set]

    parameters_list = list(stage_dict["contig"]["parameters"].keys())
    results_list += [
                     expand(output_dict["contig"] / "{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
                            genome_prefix=[config["genome_prefix"],],
                            assembly_stage=["contig",],
                            haplotype=haplotype_list,
                            parameters=parameters_list ),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["contig"],
                           haplotype=haplotype_list,
                           parameters=parameters_list),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/quast/{genome_prefix}.{assembly_stage}.{haplotype}",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["contig"],
                           haplotype=haplotype_list,
                           parameters=parameters_list),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.qv",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["contig"],
                           haplotype=haplotype_list,
                           parameters=parameters_list),
                     ] # Tested only on hifiasm


if "purge_dups" in config["stage_list"]:
    prev_stage = stage_dict["purge_dups"]["prev_stage"]
    purge_dupser_list = config["stage_coretools"]["purge_dups"]["default"]
    stage_dict["purge_dups"]["parameters"] = {}

    for purge_dupser in purge_dupser_list:
        for option_set in config["coretool_option_sets"][purge_dupser]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, purge_dupser, option_set)
                stage_dict["purge_dups"]["parameters"][parameters_label] = {}
                stage_dict["purge_dups"]["parameters"][parameters_label]["purge_dupser"] = purge_dupser
                stage_dict["purge_dups"]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][purge_dupser][option_set]

    parameters_list = list(stage_dict["purge_dups"]["parameters"].keys())
    results_list += [
                     expand(out_dir_path / "purge_dups/{parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta",
                     genome_prefix=[config["genome_prefix"], ],
                     assembly_stage=["contig"],
                     haplotype=haplotype_list,
                     parameters=parameters_list,
                     ),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}",
                        genome_prefix=[config["genome_prefix"], ],
                        assembly_stage=["purge_dups"],
                        haplotype=haplotype_list,
                        parameters=parameters_list),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/quast/{genome_prefix}.{assembly_stage}.{haplotype}",
                        genome_prefix=[config["genome_prefix"], ],
                        assembly_stage=["purge_dups"],
                        haplotype=haplotype_list,
                        parameters=parameters_list),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.qv",
                        genome_prefix=[config["genome_prefix"], ],
                        assembly_stage=["purge_dups"],
                        haplotype=haplotype_list,
                        parameters=parameters_list)
                    ]

if "hic_scaffolding" in config["stage_list"]:
    prev_stage = stage_dict["hic_scaffolding"]["prev_stage"]
    hic_scaffolder_list = config["stage_coretools"]["hic_scaffolding"]["default"]
    stage_dict["hic_scaffolding"]["parameters"] = {}

    for hic_scaffolder in hic_scaffolder_list:
        for option_set in config["coretool_option_sets"][hic_scaffolder]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, hic_scaffolder, option_set)
                stage_dict["hic_scaffolding"]["parameters"][parameters_label] = {}
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["hic_scaffolder"] = hic_scaffolder
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][hic_scaffolder][option_set]

    parameters_list = list(stage_dict["hic_scaffolding"]["parameters"].keys())

    print(stage_dict["hic_scaffolding"]["prev_stage"])

    results_list += [
                     expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{genome_prefix}.{assembly_stage}.{haplotype}.{resolution}.map.{ext}",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[prev_stage,],
                            haplotype=haplotype_list,
                            parameters=stage_dict[prev_stage]["parameters"],
                            resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                            ext=parameters["tool_options"]["pretextsnapshot"]["format"]),
                     expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{genome_prefix}.{assembly_stage}.{haplotype}.{resolution}.map.{ext}",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=["hic_scaffolding",],
                            haplotype=haplotype_list,
                            parameters=stage_dict["hic_scaffolding"]["parameters"],
                            resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                            ext=parameters["tool_options"]["pretextsnapshot"]["format"]),
                     expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=["hic_scaffolding", ],
                            haplotype=haplotype_list,
                            parameters=parameters_list),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/quast/{genome_prefix}.{assembly_stage}.{haplotype}",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["hic_scaffolding", ],
                           haplotype=haplotype_list,
                           parameters=parameters_list),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.qv",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["hic_scaffolding", ],
                           haplotype=haplotype_list,
                           parameters=parameters_list),
                    ] # TODO: implement

"""
results_dict["check_input"] = [
                               final_config_yaml,
                               final_input_yaml
                              ]

results_dict["qc"] = [*results_dict["check_input"],
                      *[expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=input_file_prefix_dict[dat_type], #list(
                                 #    map(
                                 #        lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                 #        input_filedict[dat_type])
                                 #        )
                               ) for dat_type in fastq_based_data_type_set],
                      expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                             datatype=fastq_based_data_type_set,
                             stage=["raw",]),
                      ]

results_dict["filtering"] = [*results_dict["qc"],
                             expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                                    fileprefix=input_file_prefix_dict["hifi"]) if "hifi" in fastq_based_data_type_set else [],
                             expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                                    datatype=["hifi", ],
                                    stage=["filtered", ],
                                    fileprefix=input_file_prefix_dict["hifi"],
                                    ) if "hifi" in fastq_based_data_type_set else [],
                             expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                                    datatype=["hifi"],
                                    stage=["filtered",]) if "hifi" in fastq_based_data_type_set else [], # only hifi filtration was implemented yet
                             *[[expand(output_dict["kmer"] / "{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters",
                                    datatype=[dat_type,],
                                    stage=["filtered",],
                                    kmer_tool=[kmer_tool,],
                                    kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                                    ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set]
                              ]

results_dict["contig"] = [*results_dict["filtering"],
                          #expand(output_dict["contig"] / ("{assembler}/%s.contig.{assembler}.hifi.hic.r_utg.gfa" % config["genome_name"]),
                          #       assembler=assembler_list,),
                          #output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.hifi.hic.r_utg.gfa" % config["genome_name"]),
                          expand(output_dict["contig"] / ("{assembler}/%s.contig.{assembler}.{haplotype}.fasta" % config["genome_name"]),
                                 haplotype=haplotype_list,#["p", "a"],
                                 assembler=assembler_list,),
                          expand(output_dict["assembly_qc"] /"{assembly_stage}/busco5/{assembler}/{haplotype}/",
                                 assembly_stage=["contig"],
                                 haplotype=haplotype_list,#["p", "a"],,
                                 assembler=assembler_list ,),
                          expand(output_dict["assembly_qc"] /("{assembly_stage}/quast/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}"
                                                   % config["genome_name"]),
                                 assembly_stage=["contig"],
                                 haplotype=haplotype_list,
                                 assembler=assembler_list ,),
                          expand(output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.qv" % config["genome_name"]),
                                 assembly_stage=["contig"],
                                 assembler=assembler_list ),

                          ]

results_dict["purge_dups"] = [*results_dict["contig"],
                              expand(out_dir_path / ("{assembly_stage}/{assembler}/input/%s.contig.{assembler}.{haplotype}.fasta" % config["genome_name"]),
                                     assembly_stage=["purge_dups"],
                                     assembler=assembler_list,
                                     haplotype=[primary_haplotype]),
                              expand(out_dir_path / "{assembly_stage}/{assembler}/{haplotype}/dups.bed",
                                     assembly_stage=["purge_dups"],
                                     assembler=assembler_list,
                                     haplotype=haplotype_list),
                              expand(out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"]),
                                     assembler=assembler_list,
                                     haplotype=haplotype_list,
                                     ),
                              expand(output_dict["assembly_qc"] /"{assembly_stage}/busco5/{assembler}/{haplotype}/",
                                     assembly_stage=["purge_dups"],
                                     haplotype=haplotype_list,
                                     assembler=assembler_list ,),
                              expand(output_dict["assembly_qc"] /("{assembly_stage}/quast/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}"
                                                       % config["genome_name"]),
                                     assembly_stage=["purge_dups"],
                                     haplotype=haplotype_list,
                                     assembler=assembler_list ,),
                              expand(output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.qv" % config["genome_name"]),
                                     assembly_stage=["purge_dups"],
                                     assembler=assembler_list ),
                              expand(out_dir_path  / ("purge_dups/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.fasta.bwt" % config["genome_name"]),
                                     assembly_stage=["purge_dups"],
                                     assembler=assembler_list,
                                     haplotype=haplotype_list),
                              expand(out_dir_path  / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.fai" % config["genome_name"]),
                                     assembly_stage=["purge_dups"],
                                     assembler=assembler_list,
                                     haplotype=haplotype_list),
                              expand(out_dir_path  / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.dict" % config["genome_name"]),
                                     assembly_stage=["purge_dups"],
                                     assembler=assembler_list,
                                     haplotype=haplotype_list)
                              ]
results_dict["hic_scaffolding"] = [*results_dict["purge_dups"],
                                   expand(out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{fileprefix}.bam" % config["genome_name"]),
                                          assembly_stage=["purge_dups"],
                                          haplotype=haplotype_list,
                                          assembler=assembler_list,
                                          fileprefix=input_file_prefix_dict["hic"]),
                                   expand(out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bed"  % config["genome_name"]),
                                          assembly_stage=["purge_dups"],
                                          haplotype=haplotype_list,
                                          assembler=assembler_list,),
                                   expand(out_dir_path / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.{resolution}.map.{ext}" % config["genome_name"]),
                                          assembly_stage=["purge_dups"],
                                          haplotype=haplotype_list,
                                          assembler=assembler_list,
                                          resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                          ext=parameters["tool_options"]["pretextsnapshot"]["format"]),
                                   expand(out_dir_path  / ("hic_scaffolding/{assembler}/%s.hic_scaffolding.{assembler}.{haplotype}.fasta" % config["genome_name"]),
                                          haplotype=haplotype_list,
                                          assembler=assembler_list,),
                                   expand(output_dict["assembly_qc"] /"{assembly_stage}/busco5/{assembler}/{haplotype}/",
                                          assembly_stage=["hic_scaffolding"],
                                          haplotype=haplotype_list,
                                          assembler=assembler_list ,),
                                   expand(output_dict["assembly_qc"] /("{assembly_stage}/quast/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}" % config["genome_name"]),
                                          assembly_stage=["hic_scaffolding"],
                                          haplotype=haplotype_list,
                                          assembler=assembler_list ,),
                                   expand(output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.qv" % config["genome_name"]),
                                          assembly_stage=["hic_scaffolding"],
                                          assembler=assembler_list ),
                                   expand(out_dir_path  / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.fasta.bwt" % config["genome_name"]),
                                          assembly_stage=["hic_scaffolding"],
                                          assembler=assembler_list,
                                          haplotype=haplotype_list),
                                   expand(out_dir_path  / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.fai" % config["genome_name"]),
                                          assembly_stage=["hic_scaffolding"],
                                          assembler=assembler_list,
                                          haplotype=haplotype_list),
                                   expand(out_dir_path  / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.{haplotype}.dict" % config["genome_name"]),
                                          assembly_stage=["hic_scaffolding"],
                                          assembler=assembler_list,
                                          haplotype=haplotype_list),
                                   expand(out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.{fileprefix}.bam" % config["genome_name"]),
                                          assembly_stage=["hic_scaffolding"],
                                          haplotype=haplotype_list,
                                          assembler=assembler_list,
                                          fileprefix=input_file_prefix_dict["hic"]),
                                   expand(out_dir_path / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.{resolution}.map.{ext}" % config["genome_name"]),
                                          assembly_stage=["hic_scaffolding"],
                                          haplotype=haplotype_list,
                                          assembler=assembler_list,
                                          resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                          ext=parameters["tool_options"]["pretextsnapshot"]["format"]),
                                   ]

results_dict["full"] = results_dict["hic_scaffolding"]
"""
#----

#---- Final rule ----
rule all:
    input:
        results_list
        #results_dict[config["mode"]]
#----

#---- Include section ----
include: "workflow/rules/Preprocessing/Files.smk"
include: "workflow/rules/QCFiltering/FastQC.smk"
include: "workflow/rules/QCFiltering/MultiQC.smk"
include: "workflow/rules/QCFiltering/Cutadapt.smk"
include: "workflow/rules/Kmer/Jellyfish.smk"
include: "workflow/rules/Kmer/Meryl.smk"
include: "workflow/rules/Kmer/Genomescope.smk"
include: "workflow/rules/Contigs/Hifiasm.smk"
include: "workflow/rules/Contigs/Gfatools.smk"
include: "workflow/rules/QCAssembly/BUSCO5.smk"
include: "workflow/rules/QCAssembly/Merqury.smk"
include: "workflow/rules/QCAssembly/QUAST.smk"
include: "workflow/rules/Purge_dups/Purge_dups.smk"
include: "workflow/rules/Alignment/Index.smk"
include: "workflow/rules/Alignment/Alignment.smk"
include: "workflow/rules/Alignment/Pretext.smk"
include: "workflow/rules/HiC/Salsa2.smk"

#----
