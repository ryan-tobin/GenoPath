import subprocess
import os

current_script_dir = os.path.dirname(os.path.realpath(__file__))
r_scripts_dir = os.path.join(current_script_dir)

def run_r_script(script_name):
    script_path = os.path.join(r_scripts_dir, script_name)
    print(script_path)
    command = ["Rscript", script_path]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {os.path.basename(script_path)}: {e.stderr}")
        print(f"Standard Output: {e.stdout}")


def install_QP():
    try:
        run_r_script("install_QP.R")
        print("SignatureEstimation Downloaded Successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running Comdist script: {e.output}. Execution halted.")


def install_dSig():
    try:
        run_r_script("install_dSig.R")
        print("DeconstructSigs Downloaded Successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running Comdist script: {e.output}. Execution halted.")


def install_MutPat():
    try:
        run_r_script("install_MutPat.R")
        print("MutationalPatterns Downloaded Successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running Comdist script: {e.output}. Execution halted.")
