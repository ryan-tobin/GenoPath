import argparse
import os
import shutil
import glob
import subprocess
import pandas as pd
from graphviz import Digraph
import time


def clone_git_repo(repo_url, target_dir):
    clone_dir = os.path.join(target_dir, "Meltos")
    if not os.path.isdir(clone_dir):
        subprocess.run(["git", "clone", repo_url, clone_dir], check=True)


def load_SampC(SNV):
    SNVta = pd.read_csv(SNV, sep="\t")

    print(SNVta)

    Chrom = []
    for i in SNVta["CHR"]:
        Chrom.append(i.replace("chr", ""))

    Len = len(SNVta["CHR"])

    ColLs = SNVta.columns
    SampLs = []
    for Col in ColLs:
        if Col.find(":ref") != -1:
            SampLs.append(Col.split(":")[0])

    SampC = len(SampLs) + 1

    print(SampC)

    return SampC


def order_verification(SV, OutSNV):
    print("\n-svFile (", SV, ") should have the same sample order as ", OutSNV)
    print(
        f"Checking {SV} to make sure the sample order is the same as {OutSNV}...")

    sv_df = pd.read_csv(SV, sep="\t", engine="python")
    outsnv_df = pd.read_csv(OutSNV, sep="\t", engine="python")

    sv_sample_cols = [
        col.replace("_GenomeCounts", "")
        for col in sv_df.columns
        if "_GenomeCounts" in col
    ]

    desc_index = outsnv_df.columns.tolist().index("desc")

    new_order = outsnv_df.columns[desc_index + 1:].tolist() + [
        col for col in sv_sample_cols if col + "_GenomeCounts" in outsnv_df.columns
    ]
    reordered_outsnv_df = outsnv_df.reindex(
        columns=["#chrom", "pos", "desc"] + new_order
    )

    if outsnv_df.columns.tolist() != reordered_outsnv_df.columns.tolist():
        print("Columns reordered to match SV file sample order.")
        outsnv_df = reordered_outsnv_df
    else:
        print("Columns are already in the correct order")

    sv_df.to_csv(SV, sep="\t", index=False)
    outsnv_df.to_csv(OutSNV, sep="\t", index=False)

    print(
        f"DataFrames written back to their original files: {SV} and {OutSNV}")


def run_meltos(OutSNV, SV, OutTree, SampC, target_dir):
    meltos_repo_url = "https://github.com/ih-lab/Meltos"
    clone_git_repo(meltos_repo_url, target_dir)
    print("Meltos Repository Cloned.")

    original_dir = os.getcwd()
    meltos_dir = os.path.join(target_dir, "Meltos")
    os.chdir(meltos_dir)

    # print ('java -jar Meltos.jar -treeFile '+OutTree+' -svFile '+SV+' -numSamples '+str(SampC)+' -ssnvFile '+ OutSNV)
    subprocess.run(
        [
            "java",
            "-jar",
            "Meltos.jar",
            "-treeFile",
            OutTree,
            "-svFile",
            SV,
            "-numSamples",
            str(SampC),
            "--ssnvFile",
            OutSNV
        ]
    )
