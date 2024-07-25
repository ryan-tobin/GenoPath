import os
import sys
import subprocess
import argparse
import time
import pandas as pd
import pkg_resources


def process_clone_presence(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t")

    df.iloc[:, 1:] = df.iloc[:, 1:].map(lambda x: 1 if x > 0 else 0).astype(int)

    df.to_csv(output_file, sep="\t", index=False)


def rename_hg19_to_normal(input_file, output_file):
    with open(input_file, "r") as file:
        lines = file.readlines()

    if len(lines) >= 5:
        print(f"Original line: {lines[4].strip()}")
        if "#hg19" in lines[4]:
            lines[4] = lines[4].replace("#hg19", "#Normal")
        print(f"Modified line: {lines[4].strip()}")
    else:
        print("The input file does not have enough lines.")

    with open(output_file, "w") as file:
        file.writelines(lines)


def check_networkx_version():
    required_version = "2.8"
    try:
        installed_version = pkg_resources.get_distribution("networkx").version
    except pkg_resources.DistributionNotFound:
        raise ImportError("networkx is not installed.")

    if pkg_resources.parse_version(installed_version) > pkg_resources.parse_version(
        required_version
    ):
        raise ImportError(
            f"networkx version must be <= {required_version}. Your version: {installed_version}"
        )


def delete_files(files):
    for file in files:
        try:
            os.remove(file)
            print(f"Deleted {file}")
        except OSError as e:
            print(f"Error deleting {file}: {e.strerror}")


def t_():
    time.sleep(2)


def t(n):
    time.sleep(n)


def check_input_file_columns(input_file_path):
    required_columns = ["CHR", "Position", "Wild", "Mutant", "Trinucletide"]
    column_check_passed = True
    try:
        df = pd.read_csv(input_file_path, sep='\t')
        columns = df.columns
        for col in required_columns:
            if col not in columns:
                print(f"Column '{col}' not found in the input file.")
                column_check_passed = False

        for col in columns:
            if col not in required_columns:
                if "REF" in col:
                    new_col = col.replace("REF", "ref")
                    print(f"Column '{col}' should be renamed to '{new_col}'.")
                    column_check_passed = False
                elif "ALT" in col:
                    new_col = col.replace("ALT", "alt")
                    print(f"Column '{col}' should be renamed to '{new_col}'.")
                    column_check_passed = False

    except Exception as e:
        print(f"An error occurred while checking input file columns: {e}")
        column_check_passed = False

    return column_check_passed
