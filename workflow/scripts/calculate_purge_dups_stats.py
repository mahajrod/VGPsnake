#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from pathlib import Path

from RouToolPa.Parsers.BUSCO import BUSCOtable

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--purge_dups_bed", action="store", dest="purge_dups_bed", required=True,
                    help="Bed file with purge_dups output")
parser.add_argument("-s", "--stat_file", action="store", dest="stat_file", required=True,
                    help="File wit statistics extracted from coverage file.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

artefact_set = {"HAPLOTIG", "JUNK", "REPEAT", "OVLP", "HIGHCOV"}
stat_cov_df = pd.read_csv(args.stat_file, sep="\t", header=None, names=["#scaffold", "length", "mean_cov", "median_cov"],
                          index_col=0)

purge_dups_bed_df = pd.read_csv(args.purge_dups_bed, sep="\t", header=None, names=["#scaffold", "start", "end", "type",
                                                                                   "overlapping_scaffold"], index_col=0)

purge_dups_bed_df = pd.concat([purge_dups_bed_df, stat_cov_df.loc[purge_dups_bed_df.index]], axis=1)
purge_dups_bed_df["overlap_len"] = purge_dups_bed_df["end"] - purge_dups_bed_df["start"]
purge_dups_bed_df["overlap_faction"] = purge_dups_bed_df["overlap_len"] / purge_dups_bed_df["length"]
purge_dups_bed_df.to_csv("{}.extended.bed".format(args.output_prefix), sep="\t", index=True, header=True)

stats_df = purge_dups_bed_df[["overlap_len", "type"]].groupby(by="type").agg(["count", "sum"])
print(stats_df)

stats_df.to_csv("{}.stat".format(args.output_prefix), sep="\t", index=True, header=True)

for artefact in stats_df.index.unique():
    purge_dups_bed_df.index[purge_dups_bed_df["type"] == artefact].to_series().to_csv("{0}.{1}.ids".format(args.output_prefix,
                                                                                                           str(artefact).lower()),
                                                                                      sep="\t",
                                                                                      index=False,
                                                                                      header=False)
    purge_dups_bed_df[purge_dups_bed_df["type"] == artefact].to_csv("{0}.{1}.extended.bed".format(args.output_prefix,
                                                                                                  str(artefact).lower()),
                                                                    sep="\t",
                                                                    index=False,
                                                                    header=False)

for artefact in artefact_set - set(stats_df.index.unique()):
    with open("{0}.{1}.ids".format(args.output_prefix, str(artefact).lower()), "w") as out_fd:
        pass

