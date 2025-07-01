rule get_element_counts:
    container:
        "docker://quay.io/biocontainers/mpralib:0.8.2--pyhdfd78af_0"
    conda:
        getCondaEnv("mpralib.yaml")
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[0], 50) * 1024,  # Adjust memory based on input size
    input:
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        element_counts="results/{id}/quantification/{id}.{method}.element.input.tsv.gz",
    log:
        "logs/elements/get_element_counts.{id}.{method}.log",
    benchmark:
        "benchmarks/elements/get_element_counts.{id}.{method}.tsv"
    params:
        normalize="--normalized-counts" if config["mpralib_normalized_counts"] else "",
        bc_threshold=1,
        barcodes=lambda wc: "--barcodes" if wc.method == "bcalm" else "--oligos",
    shell:
        """
        mpralib sequence-design get-counts \
        --input {input.counts} --sequence-design {input.sequence_design} \
        {params.barcodes} --all-oligos --bc-threshold {params.bc_threshold} {params.normalize} \
        --output {output.element_counts} > {log} 2>&1
        """


rule run_elements_quantification:
    container:
        "docker://visze/bcalm:latest"
    threads: 1
    resources:
        # Adjust memory based on input size
        mem_mb=lambda wc, input, attempt: calc_mem_gb(
            input[0], 450 if wc.method == "bcalm" else 50, attempt
        )
        * 1024,
    retries: 3
    input:
        element_counts="results/{id}/quantification/{id}.{method}.element.input.tsv.gz",
        labels=config.get("label_file", "UNKNOWN_LABEL_FILE"),
        script=lambda wc: getScript(f"{wc.method}_elements.R"),
    output:
        result="results/{id}/quantification/{id}.{method}.element.output.tsv.gz",
        vulcano_plot="results/{id}/quantification/{id}.{method}.element.vulcano.png",
        density_plot="results/{id}/quantification/{id}.{method}.element.density.png",
    log:
        "logs/elements/run_elements_quantification.{id}.{method}.log",
    benchmark:
        "benchmarks/elements/run_elements_quantification.{id}.{method}.tsv"
    params:
        normalize="FALSE" if config["mpralib_normalized_counts"] else "TRUE",
        control_label=config.get("control_label", "UNKNOWN_CONTROL_LABEL"),
        test_label=config.get("test_label", "UNKNOWN_TEST_LABEL"),
    shell:
        """
        Rscript {input.script} \
        --count {input.element_counts} --labels {input.labels} \
        --test-label {params.test_label} --control-label {params.control_label} \
        --normalize {params.normalize} \
        --output {output.result} \
        --output-density-plot {output.density_plot} \
        --output-vulcano-plot {output.vulcano_plot} > {log} 2>&1
        """


rule get_reporter_elements:
    container:
        "docker://quay.io/biocontainers/mpralib:0.8.2--pyhdfd78af_0"
    conda:
        getCondaEnv("mpralib.yaml")
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[0], 50) * 1024,  # Adjust memory based on input size
    input:
        quantification="results/{id}/quantification/{id}.{method}.element.output.tsv.gz",
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        "results/{id}/reporter_elements/{id}.reporter_elements.{method}.tsv.gz",
    log:
        "logs/elements/get_reporter_elements.{id}.{method}.log",
    benchmark:
        "benchmarks/elements/get_reporter_elements.{id}.{method}.tsv"
    wildcard_constraints:
        format="(reporter_elements)|(reporter_genomic_elements)",
    params:
        bc_threshold=10,
    shell:
        """
        mpralib sequence-design get-reporter-elements \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} \
        --output-reporter-elements {output} > {log} 2>&1
        """


rule get_reporter_genomic_elements:
    container:
        "docker://quay.io/biocontainers/mpralib:0.8.2--pyhdfd78af_0"
    conda:
        getCondaEnv("mpralib.yaml")
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[0], 50) * 1024,  # Adjust memory based on input size
    input:
        quantification="results/{id}/quantification/{id}.{method}.element.output.tsv.gz",
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        "results/{id}/reporter_genomic_elements/{id}.reporter_genomic_elements.{method}.bed.gz",
    log:
        "logs/elements/get_reporter_genomic_elements.{id}.{method}.log",
    benchmark:
        "benchmarks/elements/get_reporter_genomic_elements.{id}.{method}.tsv"
    params:
        bc_threshold=10,
        reference="GRCh38",
    shell:
        """
        mpralib sequence-design get-reporter-genomic-elements \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} \
        --reference {params.reference} \
        --output-reporter-genomic-elements  {output} > {log} 2>&1
        """
