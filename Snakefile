import os
import logging
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath

import yaml

import pandas as pd

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')

#---- Functions ----
def p_distance(seq_a, seq_b, seq_len):
    dist = 0
    for i in range(0, seq_len):
        if seq_a[i] != seq_b[i]:
            dist += 1
    return dist

"""
def make_fastq_lists(fastq_dir, filename_fragment_to_mark_se_reads=".se.", input_is_se=False,
                     fastq_extensions=(".fastq", ".fq")):

    fastq_dir_path = fastq_dir if isinstance(fastq_dir, PosixPath) else Path(fastq_dir)

    if not fastq_dir_path.exists():
        raise ValueError(
            "ERROR!!! Input for function 'make_fastq_lists' doesn't exist: {0}".format(str(fastq_dir_path)))

    if not fastq_dir_path.is_dir():
        raise ValueError("ERROR!!! Input for function 'make_fastq_lists' is not directory: {0}".format(str(fastq_dir_path)))

    filtered_filelist = []
    filetypes = set()
    extensions_set = set()
    for filename in fastq_dir_path.iterdir():
        suffixes = filename.suffixes()
        for extension in fastq_extensions:
            if extension in suffixes:
                break
        else:
            continue # skip file if no fastq extension was found in suffixes

        if suffixes[-1] == "gzipped":
            filetypes.add(".gz")
            extensions_set.add(suffixes[-2] + suffixes[-1])
        elif suffixes[-1] == ".bz2":
            filetypes.add("bzipped")
            extensions_set.add(suffixes[-2] + suffixes[-1])
        else:
            filetypes.add("fastq")
            extensions_set.add(suffixes[-1])

        filtered_filelist.append(filename)

    if len(filetypes) > 1:
        print("WARNING: mix of archives of different types and/or uncompressed files")

    if len(extensions_set) > 1:
        print("WARNING: mix of different extensions")

    if input_is_se:
        return filetypes, [], [], filtered_filelist

    single_end_filelist = []
    paired_end_filelist = []

    for entry in filtered_filelist:
        if filename_fragment_to_mark_se_reads in entry:
            single_end_filelist.append(entry)
        else:
            paired_end_filelist.append(entry)

    forward_filelist = paired_end_filelist[::2]
    reverse_filelist = paired_end_filelist[1:][::2]
    if len(forward_filelist) != len(reverse_filelist):
        raise ValueError(
                         "ERROR!!! Lists of forward and reverse fastqs have different length:\n"
                         "\tforward: {0}\n\treverse:{1}".format(
                                                                ",".join(list(map(str, forward_filelist))),
                                                                ",".join(list(map(str, reverse_filelist)))
                                                                )
                         )
    for forward, reverse in zip(forward_filelist, reverse_filelist):
        if len(forward.name) != len(reverse.name):
            raise ValueError("ERROR!!! Filenames of forward and reverse fastqs have different length:\n"
                             "\tforward: {0}\n\treverse:{1}".format(str(forward), str(reverse)))
        if p_distance(forward, reverse) > 1:
            raise ValueError("ERROR!!! Filenames of forward and reverse fastqs differ by more than one symbol:\n"
                             "\tforward: {0}\n\treverse:{1}".format(str(forward), str(reverse)))
    return filetypes[0], extensions_set[0], forward_filelist, reverse_filelist, single_end_filelist

"""

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
logging.info("Initialization of path variables...")
#---- Initialization of path variables for input----
input_dir_path = Path(config["input_dir"])

input_dict = {}
data_types = config["data_types"].split("_")

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
logging.info("Setting and adjusting pipeline mode...")

#pipeline_mode = config["mode"]
#starting_point = config["starting_point"]

#-------- Verification of input datatypes --------

fastq_based_data_type_set = set(data_types) & set(config["fastq_based_data"])
genome_size_estimation_data_type_set = set(config["genome_size_estimation_data"]) & fastq_based_data_type_set

logging.info("Verifying datatypes...")
for d_type in data_types:
    if d_type not in config["allowed_data_types"]:
        logging.error("Unknown data type: {0}".format(d_type))
        raise ValueError("ERROR!!! Unknown data type: {0}".format(d_type))

#--------

#----

#---- Checking input files ----
logging.info("Checking input files...")

input_filedict = {}
input_file_prefix_dict = {}

for d_type in fastq_based_data_type_set:
    input_filedict[d_type] = find_fastqs(input_dict[d_type]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict[d_type] = list(map(lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                                input_filedict[d_type]))
# check filenames of paired data
for d_type in set(config["paired_fastq_based_data"]) & fastq_based_data_type_set:
   if (len(input_filedict[d_type]) % 2) != 0:
        raise ValueError("ERROR!!! {0} fastq files seems to be unpaired or misrecognized".format(d_type))
   for forward, reverse in zip(input_filedict[d_type][::2], input_filedict[d_type][1::2]):
        print(forward, reverse)
        if p_distance(str(forward), str(reverse), len(str(forward))) > 1:
            raise ValueError("ERROR!!! Forward and reverse read files differs by more than one symbol:\n\t{0}\n\t{1}".format(str(forward),
                                                                                                                             str(reverse)))

"""
if "pacbio" in data_types:
    input_filedict["pacbio"] = find_fastqs(input_dict["pacbio"]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict["pacbio"] = list(map(lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                                input_filedict["pacbio"]))

if "nanopore" in data_types:
    input_filedict["nanopore"] = find_fastqs(input_dict["nanopore"]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict["nanopore"] =
    
if "hic" in data_types:
    input_filedict["hic"] = find_fastqs(input_dict["hic"]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict["hic"] =
    
if "lr" in data_types:
    input_filedict["lr"] = find_fastqs(input_dict["lr"]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict["lr"] =
"""

if "bionano" in data_types: # TODO: modify when input for bionano will be clear
    input_filedict["bionano"] = find_cmap(input_dict["bionano"]["dir"], cmap_extension=config["cmap_extension"])


#---- Initialize tool parameters ----
logging.info("Initializing tool parameters...")

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
    logging.info("Warning! final_kmer_counter is not in kmer_counter_list! Added...")

#check if final_kmer_length is present in parameters of final_kmer_tool
for dat_type in genome_size_estimation_data_type_set:
    if config["final_kmer_length"] not in parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"]:
        parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"].append(config["final_kmer_length"])
        logging.info("Warning! Final_kmer_length is not in parameters of final_kmer_counter! Added...")

#----

#---- Check configuration ----
#if config["mode"] in ["contig",]:
#    if len(config["kmer_counter_list"]) > 1:
#        raise ValueError("ERROR!!! Multiple kmer counter tools are not allowed in mode {0}. "
#                         "Select one.".format(config["mode"]))
#    for kmer_tool in config["kmer_counter_list"]:
#        if len(parameters["tool_options"][kmer_tool]["pacbio"]["kmer_length"]) > 1:
#            raise ValueError("ERROR!!! Multiple kmer lengths are not allowed in mode {0}. "
#                             "Select one.".format(config["mode"]))

#----

#---- Save configuration and input files ----
final_config_yaml = output_dict["config"] / "config.final.yaml"
final_input_yaml = output_dict["config"] / "input.final.yaml"

os.makedirs(output_dict["config"], exist_ok=True)

with open(final_config_yaml, 'w') as final_config_fd, open(final_input_yaml, 'w') as final_input_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False)
    yaml.dump(convert_posixpath2str_in_dict(input_dict), final_input_fd, default_flow_style=False)


#-------------------------------------------
localrules: all
ruleorder: create_fastq_links > fastqc

results_dict = {}

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
                             expand(output_dict["data"] / ("fastq/pacbio/filtered/{fileprefix}%s" % config["fastq_extension"]),
                                    fileprefix=input_file_prefix_dict["pacbio"]) if "pacbio" in fastq_based_data_type_set else [],
                             expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                                    datatype=["pacbio", ],
                                    stage=["filtered", ],
                                    fileprefix=input_file_prefix_dict["pacbio"],
                                    ) if "pacbio" in fastq_based_data_type_set else [],
                             expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                                    datatype=["pacbio"],
                                    stage=["filtered",]) if "pacbio" in fastq_based_data_type_set else [], # only pacbio filtration was implemented yet
                             *[[expand(output_dict["kmer"] / "{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters",
                                    datatype=[dat_type,],
                                    stage=["filtered",],
                                    kmer_tool=[kmer_tool,],
                                    kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                                    ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set]
                              ]

results_dict["contig"] = [*results_dict["filtering"],
                          output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.pacbio.hic.r_utg.gfa" % config["genome_name"]),
                          expand(output_dict["contig"] / ("{assembler}/%s.contig.{assembler}.pacbio.hic.{haplotype}_ctg.fasta" % config["genome_name"]),
                                 haplotype=["hap1.p", "hap2.p"],#["p", "a"],
                                 assembler=["hifiasm",],),
                          expand(output_dict["assembly_qc"] /("{assembly_stage}/busco5/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}"
                                                   % config["genome_name"]),
                                 assembly_stage=["contig"],
                                 haplotype=["hap1.p", "hap2.p"],#["p", "a"],,
                                 assembler=["hifiasm",],),
                          expand(output_dict["assembly_qc"] /("{assembly_stage}/quast/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic.{haplotype}"
                                                   % config["genome_name"]),
                                 assembly_stage=["contig"],
                                 haplotype=["hap1.p", "hap2.p"],
                                 assembler=["hifiasm",],),
                          expand(output_dict["assembly_qc"] /("{assembly_stage}/merqury/{assembler}/%s.{assembly_stage}.{assembler}.pacbio.hic"
                                                   % config["genome_name"]),
                                 assembly_stage=["contig"],
                                 assembler=["hifiasm",])
                          ]


#TODO: implement following modes when necessary
"""
results_dict["basecall"] =
results_dict["basecall_pacbio"] =
results_dict["basecall_hic"] =
results_dict["basecall_10x"] =
results_dict[create_map_bionano"] = 

"""
rule all:
    input:
        results_dict[config["mode"]]

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