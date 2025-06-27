import os, math

SCRIPTS_DIR = "../scripts"
ENVS_DIR = "../envs"


def getWorkflowFile(dir_name, name):
    return workflow.source_path("%s/%s" % (dir_name, name))


def getScript(name):
    return getWorkflowFile(SCRIPTS_DIR, name)


def getCondaEnv(name):
    return getWorkflowFile(ENVS_DIR, name)


def getLabelFile():
    if "label_file" in config:
        return config["label_file"]
    else:
        raise WorkflowError("You must define 'label_file' in the config.")


def calc_mem_gb(input_file, scaling_factor=10, attemt=1):
    """
    Calculate the size of the input file in MB.
    """
    file_size_mb = os.path.getsize(input_file) / (1024 * 1024)  # Convert bytes to MB
    mem_gb = math.ceil(
        (file_size_mb * scaling_factor) / 1024
    )  # Convert MB to GB and scale

    # Minimum 1 GB
    return max(mem_gb, 1) * attemt


configfile: getWorkflowFile("../../config", "config.yaml")
