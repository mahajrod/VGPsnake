from collections import OrderedDict
from copy import deepcopy
import pandas as pd
from pathlib import Path



#-- Initialization of path variables from config file --

#---- Initialization of path variables for input----


#---- Initialization of path variables for output directories ----

#------- Initialization of path variables for output variables --------




#---- Initialization path variables for resources ----

#-------- Reference path variables --------




#-------- Variant filtration path variables ---------



#-------- Annotation path variables --------

#----

#-------- Verify index availability --------

#---- Setting mode of pipeline ----
"""
pipeline_mode = config["mode"]

if pipeline_mode in ["both", "test", "index"]:
    mol_type_list = ["dna", "rna"]
elif pipeline_mode in ["rna", "index_rna"]:
    mol_type_list = ["rna"]
elif pipeline_mode in ["dna", "panel_of_normals", "index_dna", "basecall"]:
    mol_type_list = ["dna"]
else:
    raise ValueError("ERROR!!! Unknown pipeline mode! Allowed five modes: dna, rna, both, test, panel_of_normals")

working_fastq_dir_path = test_fastq_dir_path if pipeline_mode == "test" else input_fastq_dir_path # is used only if config["starting_point"] == "fastq"

if pipeline_mode == "test":
    pass
"""

# ---- Use embedded files in docker if custom are absent ----

# ----
#-------------------------------------------


 #---- Parsing library types and replicates ----


#---- Functions ----

#----
"""
localrules: all


if pipeline_mode == "index":
    rule all:
        input:
            

elif pipeline_mode == "index_dna":
    rule all:
        
elif pipeline_mode == "index_rna":
    rule all:
        

elif pipeline_mode == "basecall":
    
elif pipeline_mode == "panel_of_normals":
    rule all:
        input:
            

elif pipeline_mode in ["both", "test"]:
    rule all:
        input:
            


if pipeline_mode in ["index", "index_rna", "index_dna"]:
    include: "workflow/rules/Preprocessing/Reference.smk"
elif pipeline_mode in ["basecall", ]:
    include: "workflow/rules/BaseCall/BaseCall.smk"
else:
    include: "workflow/rules/Preprocessing/Target.smk"
    include: "workflow/rules/Preprocessing/AnnotationSources.smk"
    include: "workflow/rules/BaseCall/BaseCall.smk"
    include: "workflow/rules/QCFiltering/Cutadapt.smk"
    include: "workflow/rules/QCFiltering/FastQC.smk"
    if config["umi_type"] == "UMI_UDI":
        include: "workflow/rules/Alignment/ConsensusCall/UMI_UDI.smk"
    elif config["umi_type"] == "UMI_duplex":
        include: "workflow/rules/Alignment/ConsensusCall/UMI_duplex.smk"
    include: "workflow/rules/Alignment/ConsensusCall/Common.smk"
    include: "workflow/rules/Alignment/RawReads.smk"
    include: "workflow/rules/VariantCall/Pisces.smk"
    include: "workflow/rules/VariantCall/Mutect2.smk"
    include: "workflow/rules/VariantCall/PanelOfNormals.smk"
    include: "workflow/rules/Alignment/CountReads.smk"
    include: "workflow/rules/Annotation/Filter.smk"
    include: "workflow/rules/Annotation/SnpEff.smk"
    include:  "workflow/rules/Stats/CollectStats.smk"

"""
