SCRIPTS_DIR = "../scripts"
ENVS_DIR = "../envs"


def getWorkflowFile(dir_name, name):
    return workflow.source_path("%s/%s" % (dir_name, name))


def getScript(name):
    return getWorkflowFile(SCRIPTS_DIR, name)


def getLabelFile():
    if "label_file" in config:
        return config["label_file"]
    else:
        raise WorkflowError("You must define 'label_file' in the config.")
