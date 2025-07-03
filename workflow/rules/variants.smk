rule get_variant_counts:
    container:
        "docker://quay.io/biocontainers/mpralib:0.8.2--pyhdfd78af_0"
    conda:
        getCondaEnv("mpralib.yaml")
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[0], 75) * 1024,  # Adjust memory based on input size
    input:
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        variant_counts="results/{id}/quantification/{id}.{level}.variant.input.tsv.gz",
    log:
        "logs/variants/get_variant_counts.{id}.{level}.log",
    benchmark:
        "benchmarks/variants/get_variant_counts.{id}.{level}.tsv"
    params:
        normalize="--normalized-counts" if config["mpralib_normalized_counts"] else "",
        bc_threshold=1,
        barcodes=lambda wc: f"--{wc.level}s",
    shell:
        """
        mpralib sequence-design get-variant-counts \
        --input {input.counts} --sequence-design {input.sequence_design} \
        {params.barcodes} --bc-threshold {params.bc_threshold} {params.normalize} \
        --output {output.variant_counts} > {log} 2>&1
        """


rule get_variant_map:
    container:
        "docker://quay.io/biocontainers/mpralib:0.8.2--pyhdfd78af_0"
    conda:
        getCondaEnv("mpralib.yaml")
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[0], 10) * 1024,  # Adjust memory based on input size
    input:
        sequence_design=config["sequence_design_file"],
    output:
        variant_map="results/{id}/{id}.variant_map.tsv.gz",
    log:
        "logs/variants/get_variant_map.{id}.log",
    benchmark:
        "benchmarks/variants/get_variant_map.{id}.tsv"
    shell:
        """
        mpralib sequence-design get-variant-map \
        --sequence-design {input.sequence_design} \
        --output {output.variant_map} > {log} 2>&1
        """


rule run_variants_barcode_quantification:
    container:
        "docker://visze/bcalm:latest"
    threads: 1
    resources:
        mem_mb=lambda wc, input, attempt: calc_mem_gb(input[0], 450, attempt) * 1024,  # Adjust memory based on input size
    retries: 3
    input:
        variant_counts="results/{id}/quantification/{id}.barcode.variant.input.tsv.gz",
        variant_map="results/{id}/{id}.variant_map.tsv.gz",
        script=getScript("barcode_level_variants.R"),
    output:
        result="results/{id}/quantification/{id}.barcode.variant.output.tsv.gz",
        volcano_plot="results/{id}/quantification/{id}.barcode.variant.volcano.png",
    log:
        "logs/variants/run_variants_barcode_quantification.{id}.log",
    benchmark:
        "benchmarks/variants/run_variants_barcode_quantification.{id}.tsv"
    params:
        normalize="FALSE" if config["mpralib_normalized_counts"] else "TRUE",
    shell:
        """
        Rscript {input.script} \
        --count {input.variant_counts} --map {input.variant_map} \
        --normalize {params.normalize} \
        --output {output.result} --output-plot {output.volcano_plot} > {log} 2>&1
        """


rule run_variants_oligo_quantification:
    container:
        "docker://visze/bcalm:latest"
    threads: 1
    resources:
        mem_mb=lambda wc, input, attempt: calc_mem_gb(input[0], 70, attempt) * 1024,  # Adjust memory based on input size
    retries: 3
    input:
        variant_counts="results/{id}/quantification/{id}.oligo.variant.input.tsv.gz",
        script=getScript("oligo_level_variants.R"),
    output:
        result="results/{id}/quantification/{id}.oligo.variant.output.tsv.gz",
        vulcano_plot="results/{id}/quantification/{id}.oligo.variant.vulcano.png",
    log:
        "logs/variants/run_variants_oligo_quantification.{id}.log",
    benchmark:
        "benchmarks/variants/run_variants_oligo_quantification.{id}.tsv"
    params:
        normalize="FALSE" if config["mpralib_normalized_counts"] else "TRUE",
    shell:
        """
        Rscript {input.script} \
        --count {input.variant_counts} \
        --normalize {params.normalize} \
        --output {output.result} --output-plot {output.vulcano_plot} > {log} 2>&1
        """


rule get_reporter_variants:
    container:
        "docker://quay.io/biocontainers/mpralib:0.8.2--pyhdfd78af_0"
    conda:
        getCondaEnv("mpralib.yaml")
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[1], 75) * 1024,  # Adjust memory based on input size
    input:
        quantification="results/{id}/quantification/{id}.{level}.variant.output.tsv.gz",
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        "results/{id}/reporter_variants/{id}.reporter_variants.{level}.tsv.gz",
    wildcard_constraints:
        level="(barcode)|(oligo)",
    log:
        "logs/variants/get_reporter_variants.{id}.{level}.log",
    benchmark:
        "benchmarks/variants/get_reporter_variants.{id}.{level}.tsv"
    params:
        bc_threshold=10,
    shell:
        """
        mpralib sequence-design get-reporter-variants \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} \
        --output-reporter-variants {output} > {log} 2>&1
        """


rule get_reporter_genomic_variants:
    container:
        "docker://quay.io/biocontainers/mpralib:0.8.2--pyhdfd78af_0"
    conda:
        getCondaEnv("mpralib.yaml")
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[1], 75) * 1024,  # Adjust memory based on input size
    input:
        quantification="results/{id}/quantification/{id}.{level}.variant.output.tsv.gz",
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        "results/{id}/reporter_genomic_variants/{id}.reporter_genomic_variants.{level}.bed.gz",
    wildcard_constraints:
        level="(barcode)|(oligo)",
    log:
        "logs/variants/get_reporter_genomic_variants.{id}.{level}.log",
    benchmark:
        "benchmarks/variants/get_reporter_genomic_variants.{id}.{level}.tsv"
    params:
        bc_threshold=10,
    shell:
        """
        mpralib sequence-design get-reporter-genomic-variants \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} \
        --output-reporter-genomic-variants {output} > {log} 2>&1
        """
