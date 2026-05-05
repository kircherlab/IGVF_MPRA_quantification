import os, math

from snakemake.utils import validate

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
    mem_gb = math.ceil((file_size_mb * scaling_factor) / 1024)  # Convert MB to GB and scale

    # Minimum 1 GB
    return max(mem_gb, 1) * attemt


# Workaround: validate() is broken from Snakemake 9.5.1 to snakemake 9.14.7 in remote jobs
if version.parse(snakemake.__version__) >= version.parse("9.5.1") and version.parse(snakemake.__version__) <= version.parse(
    "9.14.7"
):
    from snakemake_interface_executor_plugins.settings import ExecMode

    # Use the global 'workflow' variable directly as recommended by Snakemake
    if workflow.remote_exec:
        old_exec_mode = workflow.exec_mode
        workflow.workflow_settings.exec_mode = ExecMode.DEFAULT
        validate(config, schema="../schemas/config.schema.yml")
        workflow.workflow_settings.exec_mode = old_exec_mode
    else:
        validate(config, schema="../schemas/config.schema.yml")
else:
    validate(config, schema="../schemas/config.schema.yml")
