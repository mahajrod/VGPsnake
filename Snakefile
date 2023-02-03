import logging
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
    return  list(fastq_dir_path.glob("*{0}".format(fastq_extension)))


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


"""
pacbio_dir_path = input_dir_path / config["pacbio_subdir"]
fastq_pacbio_dir_path = pacbio_dir_path / "fastq/"
run_pacbio_dir_path = pacbio_dir_path / "run/"

nanopore_dir_path = input_dir_path / config["nanopore_subdir"]
fastq_nanopore_dir_path = nanopore_dir_path / "fastq/"
run_nanopore_dir_path = nanopore_dir_path / "run/"

hic_dir_path = input_dir_path / config["hic_subdir"]
fastq_hic_dir_path = hic_dir_path / "fastq/"
run_hic_dir_path = hic_dir_path / "run/"

lr_dir_path = input_dir_path / config["lr_subdir"]
fastq_lr_dir_path = lr_dir_path / "fastq/"
run_lr_dir_path = lr_dir_path / "run/"


bionano_dir_path = input_dir_path / config["bionano_subdir"]
"""
#----
#---- Initialization of path variables for output ----
out_dir_path = Path(config["out_dir"])
output_dict = {}

for first_level_sub_dir in config["first_level_subdir_list"]:
    output_dict[first_level_sub_dir] = out_dir_path / config["{0}_subdir".format(first_level_sub_dir)]
"""
out_qc_dir_path = out_dir_path / config["qc_subdir"]
out_fastqc_dir_path = out_dir_path / config["qc_subdir"]

out_log_dir_path = out_dir_path / config["log_subdir"]
out_error_dir_path = out_dir_path / config["error_subdir"]
out_benchmark_dir_path = out_dir_path / config["benchmark_subdir"]
out_cluster_log_dir_path = out_dir_path / config["cluster_log_subdir"]
"""

#----
#---- Initialization path variables for resources ----
#----
#---- Setting mode of pipeline ----
logging.info("Setting and adjusting pipeline mode...")

#pipeline_mode = config["mode"]
#starting_point = config["starting_point"]



#-------- Verification of input datatypes --------

fastq_based_data_type_set = set(data_types) & set(config["fastq_based_data"])

logging.info("Verifying datatypes...")
for d_type in data_types:
    if d_type not in config["allowed_data_types"]:
        logging.error("Unknown data type: {0}".format(d_type))
        raise ValueError("ERROR!!! Unknown data type: {0}".format(d_type))

#if config["mode"] in
#--------

#----

#---- Checking input files ----
logging.info("Checking input files...")

input_filedict = {}

if "pacbio" in data_types:
    input_filedict["pacbio"] = find_fastqs(input_dict["pacbio"]["fastq_dir"], fastq_extension=config["fastq_extension"])

if "nanopore" in data_types:
    input_filedict["nanopore"] = find_fastqs(input_dict["nanopore"]["fastq_dir"], fastq_extension=config["fastq_extension"])

if "hic" in data_types:
    input_filedict["hic"] = find_fastqs(input_dict["hic"]["fastq_dir"], fastq_extension=config["fastq_extension"])

if "lr" in data_types:
   input_filedict["lr"] = find_fastqs(input_dict["lr"]["fastq_dir"], fastq_extension=config["fastq_extension"])

if "bionano" in data_types: # TODO: modify when input for bionano will be clear
    input_filedict["bionano"] = find_cmap(input_dict["bionano"]["dir"], cmap_extension=config["cmap_extension"])


#---- Initialize tool parameters ----
logging.info("Initializing tool paremeters...")

if config["parameter_set"] not in config["parameters"]:
    raise ValueError("Error!!! Unknown set of tool parameters: {0}".format(config["parameter_set"]))

copy_absent_entries(config["parameters"]["default"], config["parameters"][config["parameter_set"]]) # set default values for options absent in  "parameter_set"

for key in config["parameters"].keys(): # remove unused sets of parameters
    if key != config["parameter_set"]:
        config["parameters"].pop(key)

parameters = config["parameters"][config["parameter_set"]] # short alias for used set of parameters

#----

#---- Save configuration and input files ----
final_config_yaml = output_dict["config"] / "config.final.yaml"
final_input_yaml = output_dict["config"] / "input.final.yaml"

with open(final_config_yaml, 'w') as final_config_fd, open(final_input_yaml, 'w') as final_input_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False)
    yaml.dump(convert_posixpath2str_in_dict(input_dict), final_input_fd, default_flow_style=False)


#-------------------------------------------
localrules: all

if config["mode"] == "check_input":
    rule all:
        input:
            final_config_yaml,
            final_input_yaml

elif config["mode"] == "qc":
    rule all:
        input:
            final_config_yaml,
            final_input_yaml,
            expand(output_dict["data"] / "{datatype}/fastq/", datatype=fastq_based_data_type_set),
            *[expand(
                     output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                     datatype=[dat_type, ],
                     stage=["raw", ],
                     fileprefix=list(
                                     map(
                                         lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                         input_filedict[dat_type])
                                         )
                     ) for dat_type in fastq_based_data_type_set]

elif config["mode"] == "filtering":
    rule all:
        input:
            final_config_yaml,
            final_input_yaml

elif config["mode"] == "basecall":
    pass # TODO: implement this mode later
    #rule all:
    #    input:

elif config["mode"] == "basecall_pacbio":
    pass # TODO: implement this mode later
    #rule all:
    #    input:

elif config["mode"] == "basecall_hic":
    pass # TODO: implement this mode later
    #rule all:
    #    input:
elif config["mode"] == "basecall_10x":
    pass # TODO: implement this mode later
    #rule all:
    #    input:
elif config["mode"] == "create_map_bionano":
    pass # TODO: implement this mode later
    #rule all:
    #    input:

elif config["mode"] == "contig":
    pass # TODO: implement this mode
    #rule all:
    #    input:

include: "workflow/rules/Preprocessing/Files.smk"
include: "workflow/rules/QCFiltering/FastQC.smk"


"""
elif config["mode"] == "full":
    rule all:
        input:
            pass
"""

"""

if pipeline_mode in ["index", "index_rna", "index_dna"]:
    include: "workflow/rules/Preprocessing/Reference.smk"
elif pipeline_mode in ["basecall", ]:
    include: "workflow/rules/BaseCall/BaseCall.smk"
else:
    include: "workflow/rules/Preprocessing/Target.smk"
    include: "workflow/rules/Preprocessing/Files.smk"
    include: "workflow/rules/BaseCall/BaseCall.smk"
    include: "workflow/rules/QCFiltering/Cutadapt.smk"
    include: "workflow/rules/QCFiltering/FastQC.smk"
    if config["umi_type"] == "UMI_UDI":
        include: "workflow/rules/Alignment/ConsensusCall/UMI_UDI.smk"
    elif config["umi_type"] == "UMI_duplex":
        include: "workflow/rules/Alignment/ConsensusCall/UMI_duplex.smk"
    include: "workflow/rules/Alignment/ConsensusCall/Common.smk"
    include: "workflow/rules/Alignment/RawReads.smk"


    include: "workflow/rules/Annotation/Filter.smk"
    include: "workflow/rules/Annotation/SnpEff.smk"


"""
