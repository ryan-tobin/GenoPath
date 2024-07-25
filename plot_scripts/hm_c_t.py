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
    command = ['Rscript', script_path] + list(map(str, args))
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {os.path.basename(script_path)}: {e.stderr}")
        print(f"Standard Output: {e.stdout}")


def plot_heatmap(clone_tree, tumor_tree, presence_file, output):
    try:
        run_r_script('tree.R', clone_tree, tumor_tree, presence_file, output)
        print("Presence Heatmap with Clone and Tumor Tree plotted.")
    except subprocess.CalledProcessError as e:
        print(f"Error running Comdist script: {e.output}. Execution halted.")
