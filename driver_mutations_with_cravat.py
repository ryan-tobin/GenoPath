import subprocess
import requests
import requests
import time
import os
import zipfile
import csv
import pandas as pd
import warnings
from analysis.cancer_dict import cancer_types_dict
import glob
import argparse


def convert_to_open_cravat_format(input_file_path, output_file_path):
    proceed_with_conversion = True

    with open(input_file_path, "r") as input_file:
        input_file.readline()
        for line in input_file:
            fields = line.strip().split("\t")
            if len(fields) > 3 and fields[3] == "+":
                proceed_with_conversion = False
                break

    if proceed_with_conversion:
        with open(input_file_path, "r") as input_file, open(
            output_file_path, "w"
        ) as output_file:
            header = input_file.readline().strip().split("\t")

            chr_index = header.index("CHR") if "CHR" in header else header.index("chr")
            pos_index = header.index("POS") if "POS" in header else header.index("pos")
            ref_index = header.index("REF") if "REF" in header else header.index("ref")
            alt_index = header.index("ALT") if "ALT" in header else header.index("alt")

            line_number = 1
            for line in input_file:
                fields = line.strip().split("\t")
                chr_val = fields[chr_index]
                pos_val = fields[pos_index]
                ref_val = fields[ref_index]
                alt_val = fields[alt_index]

                formatted_line = (
                    f"{line_number}\t{chr_val}\t{pos_val}\t+\t{ref_val}\t{alt_val}\n"
                )
                output_file.write(formatted_line)
                line_number += 1
    else:
        print("Conversion skipped due to '+' found in the fourth column of a line.")


def run_cravat(input_file, cancer_type_input, output_dir):
    print("Running Open-Cravat Analysis to Find Driver Mutations...")
    time.sleep(2)
    if cancer_type_input in cancer_types_dict.keys():
        cancer_type = cancer_type_input
        chasm_option = f"chasmplus_{cancer_type}"
    elif cancer_type_input in cancer_types_dict.values():
        cancer_type = {value: key for key, value in cancer_types_dict.items()}[
            cancer_type_input
        ]
        chasm_option = f"chasmplus_{cancer_type}"
    elif cancer_type_input is None or cancer_type_input.strip() == "":
        cancer_type = "chasmplus"
        chasm_option = "chasmplus"
    else:
        print(
            f"Invalid cancer type: {cancer_type_input}. Defaulting to pan-cancer analysis."
        )
        cancer_type = "chasmplus"
        chasm_option = "chasmplus"

    print("Installing Open-Cravat Base Modules.")
    subprocess.run(["oc", "module", "install-base"], check=True)
    print("Installing Required Annotators")
    subprocess.run(
        ["oc", "module", "install", "vest", "pubmed", chasm_option], check=True
    )

    annotators_list = ["pubmed"]
    if chasm_option != "chasmplus":
        annotators_list.append(chasm_option)
    else:
        annotators_list.append(chasm_option)

    warnings.filterwarnings("ignore", category=SyntaxWarning)
    command = (
        ["oc", "run", input_file, "-l", "hg38", "-a"]
        + annotators_list
        + ["-t", "excel", "-d", output_dir]
    )

    try:
        subprocess.run(command, check=True)
        print("OpenCravat analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print("An error occurred while running OpenCravat:", e)
    except Exception as e:
        print("A general error occurred:", e)

    save_filtered_gene_pvalues_to_txt(chasm_option, input_file, output_dir)


def save_filtered_gene_pvalues_to_txt(chasm_option, input_file, output_dir):
    formatted_chasm_option = chasm_option.replace("_", " ")
    input_file_basename = os.path.splitext(os.path.basename(input_file))[0]
    excel_file_path = os.path.join(output_dir, f"{args.input_file}.xlsx")
    output_txt_file = os.path.join(
        output_dir, f"{input_file_basename}_cravat_drivers.txt"
    )

    try:
        df = pd.read_excel(excel_file_path, sheet_name="Variant", skiprows=1)
        gene_pvalue_column = "P-value"

        if gene_pvalue_column in df.columns:
            df[gene_pvalue_column] = pd.to_numeric(
                df[gene_pvalue_column], errors="coerce"
            )
            filtered_df = df[df[gene_pvalue_column] < 0.05].copy()

            if (
                "Ref Base" in filtered_df.columns
                and "Alt Base" in filtered_df.columns
                and "Chrom" in filtered_df.columns
                and "Position" in filtered_df.columns
            ):
                filtered_df["Mutation Nucleotide Change"] = (
                    filtered_df["Ref Base"].astype(str)
                    + ">"
                    + filtered_df["Alt Base"].astype(str)
                )
                filtered_df["Mutation Summary"] = (
                    filtered_df["Chrom"].astype(str)
                    + ":"
                    + filtered_df["Position"].astype(str)
                    + " "
                    + filtered_df["Mutation Nucleotide Change"]
                )

                final_df = filtered_df[
                    [
                        "Gene",
                        gene_pvalue_column,
                        "Mutation Nucleotide Change",
                        "Mutation Summary",
                    ]
                ].rename(columns={"Gene": "Driver Gene", gene_pvalue_column: "P-value"})

                with open(output_txt_file, "w") as f:
                    final_df.to_csv(f, sep="\t", index=False, header=True)
                print(f"Filtered gene p-values saved to {output_txt_file}")
            else:
                print("Required columns for processing not found.")
        else:
            print(f"Column '{gene_pvalue_column}' not found in the Excel sheet.")
    except FileNotFoundError:
        print(f"File not found: {excel_file_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Open-Cravat Analysis to Find Driver Mutations."
    )
    parser.add_argument(
        "-i", "--input_file", required=True, help="Input file path for genomic data."
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output file path for the conversion result.",
    )
    parser.add_argument(
        "-c",
        "--cancer_type",
        default=None,
        help="Cancer type input for analysis. Leave empty for pan-cancer analysis.",
    )
    args = parser.parse_args()

    input_file_path = os.path.abspath(args.input_file)
    output_file_path = os.path.abspath(args.output_file)
    output_dir = os.path.dirname(output_file_path)

    convert_to_open_cravat_format(input_file_path, output_file_path)
    run_cravat(input_file_path, args.cancer_type, output_dir)
