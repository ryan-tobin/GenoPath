import subprocess
import os
import shutil


def move_clonefinder_outputs(input_file, target_dir):
    input_file_base = os.path.splitext(os.path.basename(input_file))[0]
    output_files = [
        f"{input_file_base}snv_CloneFinder.nwk",
        f"{input_file_base}snv_summary.txt",
        f"{input_file_base}snv_CloneFinder.meg",
        f"{input_file_base}snv_CloneFinder.txt",
    ]

    input_dir = os.path.dirname(input_file)

    for output_file in output_files:
        src = os.path.join(input_dir, output_file)
        dest = os.path.join(target_dir, output_file)

        if os.path.exists(src) and not os.path.exists(dest):
            shutil.move(src, dest)
        elif os.path.exists(src) and os.path.exists(dest):
            print(f"File already exists in target directory, not moving: {output_file}")


def clone_git_repo(repo_url, target_dir):
    clone_dir = os.path.join(target_dir, "CloneFinder")
    if not os.path.isdir(clone_dir):
        subprocess.run(["git", "clone", repo_url, clone_dir], check=True)


def create_output_directory(target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)


def is_valid_line(line):
    columns = line.split("\t")
    for col in reversed(columns):
        if col.endswith(":REF") or col.endswith(":ALT"):
            continue
        else:
            break
    else:
        return False

    return True
    print("Input File is Valid")


def process_input_file(input_file, target_dir):
    with open(input_file, "r") as file:
        for line_num, line in enumerate(file, start=1):
            line = line.strip()
            if not line:
                continue

            if not is_valid_line(line):
                raise ValueError(f"Invalid format on line {line_num} in {input_file}")
    print("Input File Processed")


def run_clonefinder(mode, input_file, target_dir):
    create_output_directory(target_dir)

    clonefinder_repo_url = "https://github.com/SayakaMiura/CloneFinder"
    clone_git_repo(clonefinder_repo_url, target_dir)
    print("CloneFinder Repository Cloned.")

    original_dir = os.getcwd()
    clonefinder_dir = os.path.join(target_dir, "CloneFinder")
    os.chdir(clonefinder_dir)

    try:
        if mode == "snv":
            process_input_file(input_file, target_dir)
            subprocess.run(["python", "clonefinder.py", mode, input_file], check=True)
        else:
            raise ValueError("Unsupported mode. Supported modes: ['snv']")
    except subprocess.CalledProcessError as e:
        print(f"Error in running CloneFinder: {e}")
    finally:
        print("Moving CloneFinder Output Files to User Selected Directory")
        move_clonefinder_outputs(input_file, target_dir)
        os.chdir(original_dir)
