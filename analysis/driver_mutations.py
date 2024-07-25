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


def convert_to_cgi_format(input_file_path, output_file_path):
    with open(input_file_path, "r") as input_file, open(
        output_file_path, "w"
    ) as output_file:
        header = input_file.readline().strip().split("\t")

        chr_index = header.index(
            "CHR") if "CHR" in header else header.index("chr")
        pos_index = header.index(
            "POS") if "POS" in header else header.index("pos")
        ref_index = header.index(
            "REF") if "REF" in header else header.index("ref")
        alt_index = header.index(
            "ALT") if "ALT" in header else header.index("alt")

        output_file.write("chr\tpos\tref\talt\n")

        for line in input_file:
            fields = line.strip().split("\t")
            chr_val = fields[chr_index]
            pos_val = fields[pos_index]
            ref_val = fields[ref_index]
            alt_val = fields[alt_index]

            formatted_line = f"{chr_val}\t{pos_val}\t{ref_val}\t{alt_val}\n"
            output_file.write(formatted_line)


def submit_job(mutation_file, token, cancer_type_input, email):
    headers = {"Authorization": f"{email} {token}"}
    if cancer_type_input in cancer_types_dict.keys():
        cancer_type = cancer_type_input
    elif cancer_type_input in cancer_types_dict.values():
        cancer_type = {value: key for key, value in cancer_types_dict.items()}[
            cancer_type_input
        ]
    print(f"Cancer Type Selected (Abbreviation): {cancer_type}")
    payload = {
        "cancer_type": cancer_type,
        "title": f"{cancer_type} Run 1",
        "reference": "hg38",
    }
    with open(mutation_file, "rb") as file:
        files = {"mutations": file}
        r = requests.post(
            "https://www.cancergenomeinterpreter.org/api/v1",
            headers=headers,
            files=files,
            data=payload,
        )
        if r.status_code == 200:
            job_id = r.text.strip('"')
            print("Job submitted successfully. Job ID:", job_id)
            return job_id
        else:
            print(f"Error submitting job: {r.status_code}")
            print(r.text)
            return None


def check_job_completed(job_id, token, email):
    headers = {"Authorization": f"{email} {token}"}
    while True:
        try:
            r = requests.get(
                f"https://www.cancergenomeinterpreter.org/api/v1/{job_id}?action=logs",
                headers=headers,
                timeout=30,
            )
        except requests.exceptions.ConnectTimeout:
            print("The request timed out. Please check your connection and try again.")

        if r.status_code == 200:
            response = r.json()
            status = response.get("status")
            print(f"Current job status: {status}")
            if status == "Done":
                return True
            elif status == "Failed":
                print("Job failed to complete successfully.")
                return False
        else:
            print(f"Failed to check job status: {r.status_code}")
            return False
        time.sleep(10)


def download_results(job_id, token, email, output_file_name):
    headers = {"Authorization": f"{email} {token}"}
    r = requests.get(
        f"https://www.cancergenomeinterpreter.org/api/v1/{job_id}?action=download",
        headers=headers,
    )
    if r.status_code == 200:
        with open(output_file_name, "wb") as fd:
            fd.write(r.content)
        print("Download successful")
    else:
        print(f"Error downloading results: {r.status_code}")


def process_alterations_file(file_path, mutation_file, cancer_type, target_dir):
    cgi_sample_values = []
    print("Preview of the first few rows in the alterations.tsv file:")
    with open(file_path, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")
        row_count = 0
        for row in reader:
            if row_count < 5:
                print(row)
                row_count += 1
            if row.get("CGI-Oncogenic Prediction", "").startswith("driver"):
                gene = row["CGI-Gene"]
                protein_change = row["CGI-Protein Change"]
                ref = row["REF"]
                alt = row["ALT"]
                mutation_type = row["CGI-Type"]
                mutation_summary = row['chr'] + ":" + row['pos']
                mutation_string = (
                    f"{gene}\t{protein_change}\t{ref}>{alt}\t{mutation_type}\t{mutation_summary}"
                )
                cgi_sample_values.append(mutation_string)

        if row_count == 0:
            print(
                "The alterations.tsv file seems to be empty or not formatted as expected."
            )
        else:
            print("...")

    print(
        "\nCGI-Gene values for entries where 'CGI-Oncogeneic Prediction' starts with 'driver':"
    )
    if cgi_sample_values:
        for value in cgi_sample_values[:50]:
            print(value)
        if len(cgi_sample_values) > 50:
            print("...")
    else:
        print("No entries found where 'CGI-Protein' starts with 'driver'.")
    base_name = os.path.basename(mutation_file).replace(".txt", "")
    output_file_name = f"{base_name}_{cancer_type.replace(' ', '_')}_drivers.txt"
    output_file_path = os.path.join(target_dir, output_file_name)
    with open(output_file_path, "w") as output_file:
        output_file.write(
            "Driver Gene\tDriver Protein Change\tMutation Nucleotide Change\tMutation Type\tMutation Summary\n"
        )
        for value in cgi_sample_values:
            output_file.write(f"{value}\n")
    print(f"Driver mutations saved to {output_file_path}")


def extract_and_process_zip(zip_file_path, output_directory, mutation_file, cancer_type, target_dir):
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(output_directory)

    alterations_file_path = os.path.join(output_directory, 'alterations.tsv')
    if os.path.exists(alterations_file_path):
        process_alterations_file(
            alterations_file_path, mutation_file, cancer_type, target_dir)
    else:
        print(f"alterations.tsv file not found in {output_directory}")


def run_cgi(mutation_file, token, cancer_type_input, email, target_dir):
    print(
        "Running Cancer Genome Interpreter (CGI) Analysis to Find Driver Mutations..."
    )
    time.sleep(2)

    if cancer_type_input in cancer_types_dict.keys():
        cancer_type = cancer_type_input
    elif cancer_type_input in cancer_types_dict.values():
        cancer_type = {value: key for key, value in cancer_types_dict.items()}[
            cancer_type_input
        ]

    job_id = submit_job(mutation_file, token, cancer_type, email)
    if job_id:
        if check_job_completed(job_id, token, email):
            output_directory = os.getcwd()
            zip_file_name = f"{job_id}_CGI_results.zip"
            zip_file_path = os.path.join(output_directory, zip_file_name)
            download_results(job_id, token, email, zip_file_path)
            extract_and_process_zip(
                zip_file_path, output_directory, mutation_file, cancer_type, target_dir)
        else:
            print("Job not completed successfully.")
    else:
        print("Failed to submit job.")


def parse_csv_files(directory):
    mutation_data = {}
    for file_path in glob.glob(os.path.join(directory, "*.csv")):
        try:
            df = pd.read_csv(
                file_path, header=None, sep=",", na_values="NA", keep_default_na=False
            )
            if df.shape[1] < 3:
                print(
                    f"File {os.path.basename(file_path)} has less than 3 columns. Skipping."
                )
                continue
            for index, row in df.dropna().iterrows():
                positions = []
                raw_positions = str(row[2]).split(";")
                for raw_position in raw_positions:
                    if "-" in raw_position:
                        base, range_part = raw_position.split(":")
                        start_pos, end_pos = range_part.split("-")
                        positions.append(f"{base}:{start_pos}")
                        positions.append(f"{base}:{end_pos}")
                    else:
                        positions.append(raw_position)

                for position in positions:
                    if position:
                        file_identifier = os.path.basename(
                            file_path).split("_")[0]
                        mutation_data.setdefault(position, []).append(
                            (file_identifier, row[0], row[1])
                        )
        except Exception as e:
            print(f"Error processing file {os.path.basename(file_path)}: {e}")
    return mutation_data


def process_driver_file(file_path):
    driver_mutations = pd.read_csv(file_path, sep="\t")
    position_to_gene = {}
    for index, row in driver_mutations.iterrows():
        raw_positions = row["Mutation Summary"].split(" ")[0]
        gene_name = row["Driver Gene"]

        if "-" in raw_positions:
            start, end = raw_positions.split("-")
            base, start_pos = start.rsplit(":", 1)
            end_pos = end
            position_to_gene[f"{base}:{start_pos}"] = gene_name
            position_to_gene[f"{base}:{end_pos}"] = gene_name
        else:
            position_to_gene[raw_positions] = gene_name

    return driver_mutations, position_to_gene


def match_positions(mutation_data, driver_mutations, position_to_gene):
    matched_positions = {}
    for position, gene_name in position_to_gene.items():
        if position in mutation_data:
            matched_positions[position] = mutation_data[position]
            for file_identifier, mutation_type, count in mutation_data[position]:
                print(
                    f"Gene {gene_name} Position {position} has mutation count {count} in file {file_identifier} with mutation {mutation_type} and is a driver mutation."
                )
        else:
            print(
                f"No match found for position: {position} belonging to gene {gene_name}"
            )
    return matched_positions
