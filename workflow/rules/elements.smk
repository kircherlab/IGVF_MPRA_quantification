rule get_element_counts:
    # container:
    #   "docker://quay.io/biocontainers/mpralib:0.7.3--pyhdfd78af_0"
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
        mem_mb=lambda wc, input: calc_mem_gb(
            input[0], 450 if wc.method == "bcalm" else 50
        )
        * 1024,
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
    # container:
    #     "docker://quay.io/biocontainers/mpralib:0.7.3--pyhdfd78af_0"
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
        reporter_elements="results/{id}/{format}/{id}.{format}.{method}.tsv.gz",
    log:
        "logs/elements/get_reporter_elements.{id}.{method}.{format}.log",
    benchmark:
        "benchmarks/elements/get_reporter_elements.{id}.{method}.{format}.tsv"
    wildcard_constraints:
        format="(reporter_elements)|(reporter_genomic_elements)",
    params:
        bc_threshold=10,
        output_format=lambda wc: (
            ("get-reporter-elements", "--output-reporter-elements", "")
            if wc.format == "reporter_elements"
            else (
                "get-reporter-genomic-elements",
                "--output-reporter-genomic-elements",
                "--reference GRCh38",
            )
        ),
    shell:
        """
        mpralib sequence-design {params.output_format[0]} \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} {params.output_format[2]} \
        {params.output_format[1]}  {output.reporter_elements} > {log} 2>&1
        """
