#!/usr/bin/env python3

"""
GenoPath.py

Author: Ryan Tobin
Email: ryan.tobin@temple.edu
Created: 01-31-2024
Last Modified: 02/18/2024

Tested with Python version: 3.11.8
Tested with R version: 4.3.0
"""

import argparse
import os
import shutil
import glob
import subprocess
import pandas as pd

from analysis.helper_clonefinder import run_clonefinder, clone_git_repo
from analysis.helper_pathfinder import run_pathfinder
from analysis.helper_phylosignare import run_phylosignare
from analysis.helper_picante import run_comdist, run_comdistnt, run_unifrac
from analysis.makePhyloInput import makePhyloInput
from analysis.general_helper import (
    check_networkx_version,
    process_clone_presence,
    rename_hg19_to_normal,
    delete_files,
    t_,
    check_input_file_columns
)
from analysis.driver_mutations import (
    run_cgi,
    convert_to_cgi_format,
    parse_csv_files,
    process_driver_file,
    match_positions,
)
from plot_scripts.ps import run
from r_packages_install.install_r_packages import (
    install_dSig,
    install_MutPat,
    install_QP,
)
from analysis.helper_meltos import load_SampC, run_meltos, order_verification
from plot_scripts.meltos_plot import plot_meltos, find_most_recent_meltos_file
from plot_scripts.hm_c_t import plot_heatmap
from plot_scripts.CF_T_tree import clone_tree, tumor_tree

# exclusively used to not download requirements if the pipeline is run >1 times #
first_run_file = "first_run_done"
if not os.path.exists(first_run_file):
    subprocess.run(["pip", "install", "-r", "requirements.txt"], check=True)

    check_networkx_version()

    install_dSig()
    install_MutPat()
    install_QP()

    with open(first_run_file, "w") as f:
        f.write("First run completed.")
    print("First-time setup completed.")
else:
    print("This is not the first run. Skipping setup operations.")


def parse_initial_arguments():
    parser = argparse.ArgumentParser(
        description="Determine the process to run", add_help=False
    )
    parser.add_argument(
        "--run_process",
        nargs="+",
        default=["All"],
        help="Processes to run: CloneFinder, PathFinder, PhyloSignare, Picante, or All.",
    )
    parser.add_argument(
        "--target_dir",
        type=str,
        required=True,
        help="Target directory for output results.",
    )
    parser.add_argument(
        "--clonefinder_phylosig_input",
        type=str,
        help="Path to the CloneFinder input file, required if running PhyloSignare independently.",
    )
    parser.add_argument(
        "--pathfinder_input",
        type=str,
        help="Path to the PathFinder input file, required if running PathFinder independently",
    )
    parser.add_argument(
        "--driver_mutation_file",
        type=str,
        help="Input file for known driver mutations. Format should be branch, location, gene name",
        default=None,
    )
    parser.add_argument(
        "--ref_alt_file",
        type=str,
        help="File containing chr, pos, ref, alt data. Auto formatting will occur.",
    )
    parser.add_argument(
        "--tool",
        type=str,
        choices=["CGI"],
        help="Specify which driver mutation tool you would like to user. (CGI) ",
    )
    parser.add_argument("--email", type=str)
    parser.add_argument("--token", type=str)
    parser.add_argument("--cancer_type_input", type=str)

    args, remaining_argv = parser.parse_known_args()
    if "PhyloSignare" in args.run_process and args.driver_mutation_file is None:
        user_input = input(
            "Driver mutation file not specified. Do you want to calculate driver mutations using Cancer Genome Interpreter (CGI)? (y/n): "
        )
        if user_input.lower() == "y":
            if args.ref_alt_file is None:
                args.ref_alt_file = input(
                    "Please specify the --ref_alt_file containing chr, pos, ref, alt data: "
                )
            if args.tool is None:
                tool_choice = input(
                    "Please specify the tool to use ('cgi'): ").lower()
                while tool_choice not in ["cgi"]:
                    print("Invalid choice. Please enter 'cgi'.")
                    tool_choice = input(
                        "Please specify the tool to use ('cgi'): "
                    ).lower()
                args.tool = tool_choice

    return args, remaining_argv


def create_no_driver_file(target_dir):
    driver_file_path = os.path.join(target_dir, "no_driver_mutations.txt")
    with open(driver_file_path, "w") as f:
        f.write("Driver Gene\tP-value\tref>alt\tSequence Ontology\tMutation Summary\n")
        f.write("nodriver\t0.00\tND>ND\tnodriver\tchr0:0000000\n")
    return driver_file_path


def run_pipeline(args):
    if initial_args.tool == 'CGI' or initial_args.driver_mutation_file:
        input_file_path = args.input_file

    # ---------- CloneFinder Arguments ---------- [✓]
    clonefinder_skipped = False
    if initial_args.clonefinder_phylosig_input:
        input_file_base_name = os.path.splitext(
            os.path.basename(initial_args.clonefinder_phylosig_input)
        )[0]
    else:
        input_file_base_name = os.path.splitext(
            os.path.basename(args.input_file))[0]
    processed_clone_presence_output = os.path.join(
        target_dir, "processed_clone_presence.txt"
    )
    tree_file = os.path.join(
        target_dir, f"{input_file_base_name}snv_CloneFinder.nwk")

    files_to_check = [
        os.path.join(
            target_dir, rf"Supplementary Information\{input_file_base_name}snv_CloneFinder.nwk"),
        os.path.join(
            target_dir, rf"Supplementary Information\{input_file_base_name}snv_summary.txt"),
        os.path.join(
            target_dir, rf"Supplementary Information\{input_file_base_name}snv_CloneFinder.meg"),
        os.path.join(
            target_dir, rf"Supplementary Information\{input_file_base_name}snv_CloneFinder.txt"),
        os.path.join(
            target_dir, rf"Supplementary Information\processed_clone_presence.txt"),
    ]

    all_files_exist = all(os.path.exists(file) for file in files_to_check)

    if all_files_exist:
        print(
            "Required CloneFinder output files already exist. If you would like to run CloneFinder again, please use new data or confirm deletion of already exisiting files."
        )
        user_response = input(
            "Would you like to delete files and run a new session of CloneFinder? (Yes/No): "
        )
        if user_response.lower() == "yes":
            delete_files(files_to_check)
            clonefinder_skipped = False
        else:
            print("Skipping CloneFinder Analysis...")
            clonefinder_skipped = True
    else:
        clonefinder_skipped = False

    if not clonefinder_skipped:
        if "All" in processes or "CloneFinder" in processes:
            print("\n----- Running CloneFinder Analysis -----\n")
            input_file_base_name = os.path.splitext(os.path.basename(args.input_file))[
                0
            ]
            print(rf"The input file name is: {input_file_base_name}", "\n")
            t_()
            run_clonefinder(args.mode, args.input_file, args.target_dir)
            print("Processing Output Files for Use in PathFinder...")
            clonefinder_output_file = (
                rf"{args.target_dir}\{input_file_base_name}snv_CloneFinder.meg"
            )
            print("Renaming hg19 to Normal...")
            t_()
            processed_clonefinder_output_file = (
                rf"{args.target_dir}\{input_file_base_name}snv_PathFinder_processed.meg"
            )
            rename_hg19_to_normal(
                clonefinder_output_file, processed_clonefinder_output_file
            )

            print("Clone .meg File Processed...")
            clonefinder_clone_presence_file = (
                rf"{args.target_dir}\{input_file_base_name}snv_CloneFinder.txt"
            )

            processed_clone_presence_output = (
                rf"{args.target_dir}\processed_clone_presence.txt"
            )
            process_clone_presence(
                clonefinder_clone_presence_file, processed_clone_presence_output
            )
            print("Clone Presence File Processeed")

            aln = processed_clonefinder_output_file
            if "All" in processes:
                print("CloneFinder Analysis Complete. Moving on...")
            elif "CloneFinder" in processes:
                print("CloneFinder Analysis Complete. Finishing up and exiting...")
            elif "PathFinder" in processes and "CloneFinder" in processes:
                print(
                    "CloneFinder Analysis Complete. Moving onto PathFinder Analysis..."
                )
            t_()
    else:
        if "PathFinder" in processes:
            print("Moving onto PathFinder Analysis...")
        elif "PhyloSignare" in processes:
            print("Moving onto PhyloSignare Analysis")
        elif "Picante" in processes:
            print("Moving onto Picante Analysis...")
        t_()

    # ---------- PathFinder Arguments ---------- [✓]
    if "All" in processes or "PathFinder" in processes:
        print("\n----- Running PathFinder Analysis -----\n")
        aln = rf"{args.target_dir}\{input_file_base_name}snv_PathFinder_processed.meg"
        run_pathfinder(
            aln,
            processed_clone_presence_output,
            args.primary,
            args.max_graphs_per_tree,
            args.target_dir,
        )

        scratch_dir_pattern = os.path.join(args.target_dir, "scratch*")

        scratch_dirs = glob.glob(scratch_dir_pattern)

        scratch_dirs_sorted = sorted(
            scratch_dirs, key=os.path.getmtime, reverse=True)

        pathfinder_results_dir = os.path.join(
            args.target_dir, "PathFinder_Results")

        if not os.path.exists(pathfinder_results_dir):
            if scratch_dirs_sorted:
                scratch_dir = scratch_dirs_sorted[0]
                print(
                    f"Using the most recently created scratch directory: {scratch_dir}"
                )
                os.rename(scratch_dir, pathfinder_results_dir)
                print(
                    f"'scratch' folder renamed to 'PathFinder_Results' in {args.target_dir}"
                )
            else:
                print(
                    "No scratch directory found, and no existing 'PathFinder_Results' folder."
                )
        else:
            print(
                f"Using the existing 'PathFinder_Results' directory in {args.target_dir}"
            )

        t_()

        if "All" in processes:
            print("PathFinder Analysis Complete. Moving on...")
        elif "PathFinder" in processes:
            print("PathFinder Analysis Complete. Finishing up and exiting...")
        elif "PhyloSignare" in processes and "PathFinder" in processes:
            print("PathFinder Analysis Complete. Moving onto PhyloSignare Analysis...")
        t_()

    # ---------- PhyloSignare Arguments ---------- [✓]
    if "All" in processes or "PhyloSignare" in processes:
        print("\n----- Creating PhyloSignare Input File -----\n")
        print("Creating PhyloSignare Input...")
        t_()
        if "All" in processes or 'PhyloSignare' in processes:
            aln_phylosig = args.input_file
        else:
            aln_phylosig = initial_args.clonefinder_phylosig_input
        makePhyloInput(aln_phylosig, target_dir)
        input_file_base_name = os.path.splitext(
            os.path.basename(aln_phylosig))[0]
        original_output_file_path = os.path.join(
            os.path.dirname(aln_phylosig), f"{input_file_base_name}_PSF.input"
        )
        target_file_path = os.path.join(
            target_dir, f"{input_file_base_name}_PSF.input")

        try:
            shutil.move(original_output_file_path, target_file_path)
            print(f"PhyloSignare Input File moved to: {target_file_path}")
        except Exception as e:
            print(f"Error moving PhyloSignare input file: {e}")

        original_folder_path = os.path.join(
            os.path.dirname(aln_phylosig), input_file_base_name
        )
        target_folder_path = os.path.join(target_dir, input_file_base_name)

        try:
            if os.path.exists(original_folder_path) and os.path.isdir(
                original_folder_path
            ):
                shutil.move(original_folder_path, target_folder_path)
                print(
                    f"Folder {original_folder_path} moved to: {target_folder_path}")
                print(
                    "PhyloSignare Input File Created and Mutation Count Folder Moved. Moving onto analysis..."
                )
            else:
                print(f"No folder found at {original_folder_path} to move.")
        except Exception as e:
            print(f"Error moving folder: {e}")

        print("\n----- Running PhyloSignare Analysis -----\n")
        run_phylosignare(target_file_path, args.control_file, target_dir)
        # ---------- Post-Process Plotting ----------
        target_file_base_name = os.path.splitext(
            os.path.basename(target_file_path))[0]
        signature_files_dir = None

        phylosignare_dirs = glob.glob(
            os.path.join(target_dir, "*-PhyloSignare"))

        for dir_path in phylosignare_dirs:
            potential_file = os.path.join(dir_path, "PhyloSignare.txt")
            if os.path.isfile(potential_file):
                signature_files_dir = potential_file
                break

        if signature_files_dir:
            try:
                df = pd.read_csv(signature_files_dir, sep="\t")
                print(f"Successfully read {signature_files_dir}")
            except Exception as e:
                print(f"Error reading {signature_files_dir}: {e}")
        else:
            print("PhyloSignare.txt not found in any -PhyloSignare directory.")
        csv_files_dir = os.path.join(f"{target_dir}/{input_file_base_name}")
        if phylosignare_dirs:
            first_dir_path = phylosignare_dirs[0]
            summary_file_path = os.path.join(first_dir_path, "Summary.txt")
        else:
            print("No PhyloSignare directories found.")
        bar_plots_output_dir = os.path.join(
            target_dir, "combined_bar_plots.png")
        tree_output = os.path.join(
            target_dir, "Phylogenetic_Tree_complete.png")
        tree_output_driver = os.path.join(
            target_dir, "Phylogenetic_Tree_driver.png")
        if initial_args.tool == "CGI":
            convert_to_cgi_format(
                initial_args.ref_alt_file,
                os.path.join(
                    target_dir, f"{input_file_base_name}_cgi_mutation_file.txt"
                ),
            )
            mutation_file = os.path.join(
                target_dir, f"{input_file_base_name}_cgi_mutation_file.txt"
            )
            run_cgi(
                mutation_file,
                initial_args.token,
                initial_args.cancer_type_input,
                initial_args.email,
                target_dir
            )
            driver_file = os.path.join(
                target_dir,
                f"{input_file_base_name}_cgi_mutation_file_{initial_args.cancer_type_input.replace(' ', '_')}_drivers.txt",
            )
            mutation_counts = parse_csv_files(csv_files_dir)

            driver_mutations, position_to_gene = process_driver_file(
                driver_file)

            position_to_gene_name = {}
            for index, row in driver_mutations.iterrows():
                mutation_summary = row["Mutation Summary"]
                base_position = mutation_summary.split()[0]
                if "-" in base_position:
                    base_position = base_position.split("-")[0]
                position_to_gene_name[base_position] = row["Driver Gene"]

            matched_positions = match_positions(
                mutation_counts, driver_mutations, position_to_gene
            )
            print("Matched positions:")
            for position_key, matches in matched_positions.items():
                print(f"{position_key}: {matches}")

            run(
                target_dir,
                csv_files_dir,
                signature_files_dir,
                args.control_file,
                summary_file_path,
                input_file_base_name,
                tree_output,
                tree_output_driver,
                bar_plots_output_dir,
                matched_positions,
                driver_mutations,
                position_to_gene,
                position_to_gene_name,
            )
        elif initial_args.driver_mutation_file is not None:
            driver_file = initial_args.driver_mutation_file
            mutation_counts = parse_csv_files(csv_files_dir)

            driver_mutations, position_to_gene = process_driver_file(
                driver_file)

            position_to_gene_name = {}
            for index, row in driver_mutations.iterrows():
                mutation_summary = row["Mutation Summary"]
                base_position = mutation_summary.split()[0]
                if "-" in base_position:
                    base_position = base_position.split("-")[0]
                position_to_gene_name[base_position] = row["Driver Gene"]

            matched_positions = match_positions(
                mutation_counts, driver_mutations, position_to_gene
            )
            print("Matched positions:")
            for position_key, matches in matched_positions.items():
                print(f"{position_key}: {matches}")

            run(
                target_dir,
                csv_files_dir,
                signature_files_dir,
                args.control_file,
                summary_file_path,
                input_file_base_name,
                tree_output,
                tree_output_driver,
                bar_plots_output_dir,
                matched_positions,
                driver_mutations,
                position_to_gene,
                position_to_gene_name,
            )

        elif initial_args.driver_mutation_file is None:
            driver_file = create_no_driver_file(target_dir)
            mutation_counts = parse_csv_files(csv_files_dir)

            driver_mutations, position_to_gene = process_driver_file(
                driver_file)

            position_to_gene_name = {}
            for index, row in driver_mutations.iterrows():
                mutation_summary = row["Mutation Summary"]
                base_position = mutation_summary.split()[0]
                if "-" in base_position:
                    base_position = base_position.split("-")[0]
                position_to_gene_name[base_position] = row["Driver Gene"]

            matched_positions = match_positions(
                mutation_counts, driver_mutations, position_to_gene
            )
            print("Matched positions:")
            for position_key, matches in matched_positions.items():
                print(f"{position_key}: {matches}")

            run(
                target_dir,
                csv_files_dir,
                signature_files_dir,
                args.control_file,
                summary_file_path,
                input_file_base_name,
                tree_output,
                tree_output_driver,
                bar_plots_output_dir,
                matched_positions,
                driver_mutations,
                position_to_gene,
                position_to_gene_name,
            )

        if "All" in processes:
            print("PhyloSignare Analysis Complete. Graphing is done...")
        elif "PhyloSignare" in processes:
            print("PhyloSignare Analysis Complete. Graphing is done...")
        elif "Picante" in processes and "PhyloSignare" in processes:
            print("PhyloSignare Analysis Complete. Graphing is done...")
        t_()

    if 'All' in processes or 'Picante_all' in processes or 'Picante_comdist' in processes or 'Picante_comdistnt' in processes or 'Picante_Unifrac' in processes:
        if args.abundance_weighted == "TRUE":
            COMDIST_OUTPUT_PDF = os.path.join(
                target_dir, f"{input_file_base_name}_comdist_weighted.pdf"
            )
            COMDISTNT_OUTPUT_PDF = os.path.join(
                target_dir, f"{input_file_base_name}_comdistnt_weighted.pdf"
            )
            UNIFRAC_OUTPUT_PDF = os.path.join(
                target_dir, f"{input_file_base_name}_unifrac_weighted.pdf"
            )
        else:
            COMDIST_OUTPUT_PDF = os.path.join(
                target_dir, f"{input_file_base_name}_comdist.pdf"
            )
            COMDISTNT_OUTPUT_PDF = os.path.join(
                target_dir, f"{input_file_base_name}_comdistnt.pdf"
            )
            UNIFRAC_OUTPUT_PDF = os.path.join(
                target_dir, f"{input_file_base_name}_unifrac.pdf"
            )

    if clonefinder_skipped:
        comm = processed_clone_presence_output
        tree = tree_file
    else:
        comm = processed_clone_presence_output
        tree = tree_file

    COMDIST_NEWICK_OUT = os.path.join(
        target_dir, f"{input_file_base_name}_comdist_newick.txt"
    )
    COMDISTNT_NEWICK_OUT = os.path.join(
        target_dir, f"{input_file_base_name}_comdistnt_newick.txt"
    )
    UNIFRAC_NEWICK_OUT = os.path.join(
        target_dir, f"{input_file_base_name}_unifrac_newick.txt"
    )

    # ---------- Picante Comdist Arguments ---------- [✓]
    if "Picante_comdist" in processes:
        print("\n----- Running Picante Analysis -----\n")
        print(comm)
        print(tree)
        print(COMDISTNT_OUTPUT_PDF)
        run_comdist(comm, tree, args.abundance_weighted,
                    COMDIST_OUTPUT_PDF, COMDIST_NEWICK_OUT)
        if "All" or "Picante" in processes:
            print("Picante Comdist Analysis Complete. Finishing up and exiting...")
        t_()

    # ---------- Picante Comdistnt Arguments ---------- [✓]
    if "Picante_comdistnt" in processes:
        print("\n----- Running Picante Comdistnt Analysis -----\n")
        run_comdistnt(comm, tree, args.abundance_weighted,
                      COMDISTNT_OUTPUT_PDF, COMDISTNT_NEWICK_OUT)
        if "All" or "Picante_comdisnt" in processes:
            print("Picante Comdistnt Analysis Complete. Finishing up and exiting...")
        t_()

    # ---------- Picante Unifrac Arguments ---------- [✓]
    if "Picante_Unifrac" in processes:
        print("\n----- Running Picante Unifrac Analysis -----\n")
        run_unifrac(comm, tree, UNIFRAC_OUTPUT_PDF, UNIFRAC_NEWICK_OUT)
        if "All" or "Picante_unifrac" in processes:
            print("Picante Unifrac Analysis Complete. Finishing up and exiting...")
        t_()

    # ---------- Picante Arguments ---------- [✓]
    if "All" in processes or "Picante_all" in processes:
        print("\n----- Running Picante Analysis (Comdist, Comdisnt, Unifrac) -----\n")
        run_comdist(comm, tree, args.abundance_weighted,
                    COMDIST_OUTPUT_PDF, COMDIST_NEWICK_OUT)
        t_()
        run_comdistnt(comm, tree, args.abundance_weighted,
                      COMDISTNT_OUTPUT_PDF, COMDISTNT_NEWICK_OUT)
        t_()
        run_unifrac(comm, tree, UNIFRAC_OUTPUT_PDF, UNIFRAC_NEWICK_OUT)
        t_()
        if "All" or "Picante_all" in processes:
            print(
                "Picante Analysis (Comdist, Comdisnt, Unifrac) Complete. Finishing up and exiting..."
            )
            t_()

    # ---------- Meltos ----------
    if "Meltos" in processes:
        print("\n----- Creating Meltos Input File ------\n")
        if "CloneFinder" in processes:
            SNV = args.input_file
        if "PathFinder" in processes and "CloneFinder" not in processes:
            SNV = initial_args.pathfinder_input
        if (
            "PhyloSignare" in processes
            and "CloneFinder" not in processes
            and "PathFinder" not in processes
        ):
            SNV = initial_args.clonefinder_phylosig_input

        SV = args.sv_file
        input_file_base_name = os.path.splitext(os.path.basename(SNV))[0]
        print(f"Input file for base name: {input_file_base_name}")

        print("Creating Meltos In File...")
        subprocess.run(
            ["python", "CloneFinder2MeltosIn.py", SNV], check=True, capture_output=True
        )
        print("Meltos input file created.")

        OutSNV = os.path.join(
            target_dir, f"{input_file_base_name}snv_CloneFinderSNV.txt")
        print(f"SNV file to be used for Meltos is: {OutSNV}")

        OutTree = os.path.join(
            target_dir, f"{input_file_base_name}snv_CloneFinderTree.txt")

        SampC = load_SampC(SNV)
        order_verification(SV, OutSNV)

        run_meltos(OutSNV, SV, OutTree, SampC, target_dir)

        meltos_output_path = find_most_recent_meltos_file(target_dir)
        node_clone_mapping_file = os.path.join(
            target_dir, f"{input_file_base_name}snv_CloneFinderCloneID.txt"
        )

        print(f"Most recent .meltos file: {meltos_output_path}")
        print(f"Node Clone Mapping File loaded as: {node_clone_mapping_file}")

        if initial_args.tool == "CGI":
            driver_file_path = os.path.join(
                target_dir,
                f"{input_file_base_name}_cgi_mutation_file_{initial_args.cancer_type_input.replace(' ', '_')}_drivers.txt",
            )

        elif initial_args.driver_mutation_file is None:
            driver_file_path = create_no_driver_file(target_dir)

        else:
            driver_file_path = initial_args.driver_mutation_file

        plot_meltos(target_dir, SV, driver_file_path, input_file_base_name)

    # ---------- Heatmap ----------
    if 'CloneFinder' in processes and 'Picante_all' in processes or 'Picante_comdist' in processes or 'Picante_comdistnt' in processes or 'Picante_Unifrac' in processes:
        input_file_base_name = os.path.splitext(
            os.path.basename(args.input_file))[0]
        ct = os.path.join(
            target_dir, f"{input_file_base_name}snv_CloneFinder.nwk")
        presence_file = os.path.join(
            target_dir, f"{input_file_base_name}snv_CloneFinder.txt")
        presence = pd.read_csv(presence_file, sep='\t', index_col=0)
        output = os.path.join(target_dir, "clone_tumor_heatmap.png")

        if "Normal" not in presence.columns:
            presence['Normal'] = 0.0

        presence.to_csv(presence_file, sep='\t')

    if 'All' in processes:
        print("\n----- Plotting Heatmap ------\n")
        ct = os.path.join(
            target_dir, f"{input_file_base_name}snv_CloneFinder.nwk")
        presence_file = os.path.join(
            target_dir, f"{input_file_base_name}snv_CloneFinder.txt")
        presence = pd.read_csv(presence_file, sep='\t', index_col=0)
        output = os.path.join(target_dir, "clone_tumor_heatmap.png")

        if "Normal" not in presence.columns:
            presence['Normal'] = 0.0
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdistnt_newick.txt")

        plot_heatmap(ct, tt, presence_file, output)
        print("Heatmap plotted successfully...")

    if "CloneFinder" in processes and 'Picante_all' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdistnt_newick.txt")

        plot_heatmap(ct, tt, presence_file, output)
        print("Heatmap plotted successfully...")

    if "CloneFinder" in processes and 'Picante_comdist' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdist_newick.txt")

        plot_heatmap(ct, tt, presence_file, output)
        print("Heatmap plotted successfully...")

    if "CloneFinder" in processes and 'Picante_comdistnt' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdistnt_newick.txt")

        plot_heatmap(ct, tt, presence_file, output)
        print("Heatmap plotted successfully...")

    if "CloneFinder" in processes and 'Picante_unifrac' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_unifrac_newick.txt")

        plot_heatmap(ct, tt, presence_file, output)
        print("Heatmap plotted successfully...")

    # ---------- Clone, Tumor Tree vs Presence ----------
    clone_path = os.path.join(
        target_dir, f"{input_file_base_name}_clonevpresence.png")
    tumor_path = os.path.join(
        target_dir, f"{input_file_base_name}_treevpresence.png")
    if 'All' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdistnt_newick.txt")

        clone_tree(presence_file, ct, clone_path)
        tumor_tree(presence_file, tt, tumor_path)
        print("Clone, Tumor Phylo Tree vs Presence Matrix Plotted Successfully...")

    if "CloneFinder" in processes and 'Picante_all' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdistnt_newick.txt")

        clone_tree(presence_file, ct, clone_path)
        tumor_tree(presence_file, tt, tumor_path)
        print("Clone, Tumor Phylo Tree vs Presence Matrix Plotted Successfully...")

    if "CloneFinder" in processes and 'Picante_comdist' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdist_newick.txt")

        clone_tree(presence_file, ct, clone_path)
        tumor_tree(presence_file, tt, tumor_path)
        print("Clone, Tumor Phylo Tree vs Presence Matrix Plotted Successfully...")

    if "CloneFinder" in processes and 'Picante_comdistnt' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_comdistnt_newick.txt")

        clone_tree(presence_file, ct, clone_path)
        tumor_tree(presence_file, tt, tumor_path)
        print("Clone, Tumor Phylo Tree vs Presence Matrix Plotted Successfully...")

    if "CloneFinder" in processes and 'Picante_unifrac' in processes:
        print("\n----- Plotting Heatmap ------\n")
        tt = os.path.join(
            target_dir, f"{input_file_base_name}_unifrac_newick.txt")

        clone_tree(presence_file, ct, clone_path)
        tumor_tree(presence_file, tt, tumor_path)
        print("Clone, Tumor Phylo Tree vs Presence Matrix Plotted Successfully...")


# ---------- Full Arguments ----------
def parse_full_arguments(processes, remaining_argv):
    parser = argparse.ArgumentParser(description="Run the analysis pipeline.")
    if "All" in processes or "CloneFinder" in processes:
        parser.add_argument(
            "mode",
            type=str,
            help='Operation mode for CloneFinder. **IMPORTANT**: Should always be "snv"',
        )
        parser.add_argument(
            "input_file", type=str, help="Path to the main input file for CloneFinder."
        )

    if "All" in processes or "PathFinder" in processes:
        parser.add_argument(
            "--primary",
            type=str,
            help="Optional argument if primary tumor is known for PathFinder",
        )
        parser.add_argument(
            "--max_graphs_per_tree",
            type=int,
            help="Optional argument to lower computation t_ for PathFinder, recommended to be 30.",
        )

    if "All" in processes or "PhyloSignare" in processes:
        parser.add_argument(
            "control_file", type=str, help="Control File for PhyloSignare Input"
        )

    if (
        "All" in processes
        or "Picante_comdist" in processes
        or "Picante_comdistnt" in processes
        or "Picante_all" in processes
    ):
        parser.add_argument(
            "--abundance_weighted",
            default="TRUE",
            type=str,
            help="This value by default is set to TRUE. If you would like to turn off abundance weighting, set to FALSE",
        )

    if "Meltos" in processes:
        parser.add_argument("--sv_file", type=str,
                            help="Path to SV file for Meltos")

    args = parser.parse_args(remaining_argv)
    return args


if __name__ == "__main__":
    initial_args, remaining_argv = parse_initial_arguments()
    processes = initial_args.run_process
    target_dir = initial_args.target_dir
    aln_phylosig = initial_args.clonefinder_phylosig_input

    args = parse_full_arguments(processes, remaining_argv)

    args.target_dir = target_dir

    run_pipeline(args)
