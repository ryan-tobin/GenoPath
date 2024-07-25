import subprocess
import os

current_script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_script_dir)
r_scripts_dir = os.path.join(parent_dir, "r_scripts")


def run_r_script(script_name, *args):
    print(f"Current directory: {current_script_dir}")
    print(f"Parent directory: {parent_dir}")
    print(f"R Scripts directory: {r_scripts_dir}")
    script_path = os.path.join(r_scripts_dir, script_name)
    command = ["Rscript", script_path] + list(map(str, args))
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {os.path.basename(script_path)}: {e.stderr}")
        print(f"Standard Output: {e.stdout}")


def run_comdist(comm, tree, abundance_weighted, output, newick_output):
    try:
        run_r_script("comdist.R", comm, tree,
                     str(abundance_weighted), output, newick_output)
        print("Comdist Analysis Complete")
    except subprocess.CalledProcessError as e:
        print(f"Error running Comdist script: {e.output}. Execution halted.")


def run_comdistnt(comm, tree, abundance_weighted, output, newick_output):
    try:
        run_r_script("comdistnt.R", comm, tree,
                     str(abundance_weighted), output, newick_output)
        print("Comdistnt Analysis Complete")
    except subprocess.CalledProcessError as e:
        print(f"Error running Comdistnt script: {e.output}. Execution halted.")


def run_unifrac(comm, tree, output, newick_output):
    try:
        run_r_script("unifrac.R", comm, tree, output, newick_output)
        print("Unifrac Analysis Complete")
    except subprocess.CalledProcessError as e:
        print(f"Error running UniFrac script: {e.output}. Execution halted.")
