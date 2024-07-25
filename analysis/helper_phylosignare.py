import subprocess
import os


def clone_git_repo(repo_url, target_dir):
    clone_dir = os.path.join(target_dir, "PhyloSignare")
    if not os.path.exists(clone_dir):
        subprocess.run(["git", "clone", repo_url, clone_dir], check=True)


def create_output_directory(target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)


def run_phylosignare(target_file_path, control_file, target_dir):
    create_output_directory(target_dir)

    phylosignare_repo_url = "https://github.com/SayakaMiura/PhyloSignare"
    clone_git_repo(phylosignare_repo_url, target_dir)
    print("Repository Clone.")

    phylosignare_dir = os.path.join(target_dir, "PhyloSignare")
    original_dir = os.getcwd()
    os.chdir(phylosignare_dir)

    try:
        subprocess.run(
            ["python", "phylosignare.py", target_file_path, control_file], check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error in running PhyloSignare: {e}")
    finally:
        os.chdir(original_dir)
