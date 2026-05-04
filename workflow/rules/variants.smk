rule get_variant_counts:
    input:
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        variant_counts="results/{id}/quantification/{id}.{level}.variant.input.tsv.gz",
    log:
        "logs/variants/get_variant_counts.{id}.{level}.log",
    benchmark:
        "benchmarks/variants/get_variant_counts.{id}.{level}.tsv"
    conda:
        getCondaEnv("mpralib.yaml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.4--pyhdfd78af_0"
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[0], 75) * 1024,  # Adjust memory based on input size
    params:
        normalize="--normalized-counts" if config["mpralib_normalized_counts"] else "",
        bc_threshold=1,
        barcodes=lambda wc: f"--{wc.level}s",
        scaling_factor=config.get("scaling_factor", 1e9),
        pseudo_count=config.get("pseudo_count", 1),
    shell:
        """
        mpralib combine get-variant-counts \
        --input {input.counts} --sequence-design {input.sequence_design} \
        {params.barcodes} --bc-threshold {params.bc_threshold} {params.normalize} \
        --scaling-factor {params.scaling_factor} \
        --pseudo-count {params.pseudo_count} \
        --output {output.variant_counts} > {log} 2>&1
        """


rule get_variant_map:
    input:
        sequence_design=config["sequence_design_file"],
    output:
        variant_map="results/{id}/{id}.variant_map.tsv.gz",
    log:
        "logs/variants/get_variant_map.{id}.log",
    benchmark:
        "benchmarks/variants/get_variant_map.{id}.tsv"
    conda:
        getCondaEnv("mpralib.yaml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.4--pyhdfd78af_0"
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[0], 10) * 1024,  # Adjust memory based on input size
    shell:
        """
        mpralib combine get-variant-map \
        --sequence-design {input.sequence_design} \
        --output {output.variant_map} > {log} 2>&1
        """


rule run_variants_barcode_quantification:
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
    retries: 3
    conda:
        getCondaEnv("bcalm.yaml")
    container:
        "docker://visze/bcalm:latest"
    threads: 1
    resources:
        mem_mb=lambda wc, input, attempt: calc_mem_gb(input[0], 450, attempt) * 1024,  # Adjust memory based on input size
    params:
        normalize="FALSE" if config["mpralib_normalized_counts"] else "TRUE",
        normalize_size=config.get("scaling_factor", 1e9),
    shell:
        """
        Rscript {input.script} \
        --count {input.variant_counts} --map {input.variant_map} \
        --normalize {params.normalize} --normalize-size {params.normalize_size} \
        --output {output.result} --output-plot {output.volcano_plot} > {log} 2>&1
        """


rule run_variants_oligo_quantification:
    input:
        variant_counts="results/{id}/quantification/{id}.oligo.variant.input.tsv.gz",
        script=getScript("oligo_level_variants.R"),
    output:
        result="results/{id}/quantification/{id}.oligo.variant.output.tsv.gz",
        volcano_plot="results/{id}/quantification/{id}.oligo.variant.volcano.png",
    log:
        "logs/variants/run_variants_oligo_quantification.{id}.log",
    benchmark:
        "benchmarks/variants/run_variants_oligo_quantification.{id}.tsv"
    retries: 3
    conda:
        getCondaEnv("bcalm.yaml")
    container:
        "docker://visze/bcalm:latest"
    threads: 1
    resources:
        mem_mb=lambda wc, input, attempt: calc_mem_gb(input[0], 70, attempt) * 1024,  # Adjust memory based on input size
    params:
        normalize="FALSE" if config["mpralib_normalized_counts"] else "TRUE",
        normalize_size=config.get("scaling_factor", 1e9),
    shell:
        """
        Rscript {input.script} \
        --count {input.variant_counts} \
        --normalize {params.normalize} --normalize-size {params.normalize_size} \
        --output {output.result} --output-plot {output.volcano_plot} > {log} 2>&1
        """


rule get_reporter_variants:
    input:
        quantification="results/{id}/quantification/{id}.{level}.variant.output.tsv.gz",
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        "results/{id}/reporter_variants/{id}.reporter_variants.{level}.tsv.gz",
    log:
        "logs/variants/get_reporter_variants.{id}.{level}.log",
    benchmark:
        "benchmarks/variants/get_reporter_variants.{id}.{level}.tsv"
    wildcard_constraints:
        level="(barcode)|(oligo)",
    conda:
        getCondaEnv("mpralib.yaml")
    container:
        "docker://quay.io/biocontainers/mpralib:0.10.4--pyhdfd78af_0"
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[1], 75) * 1024,  # Adjust memory based on input size
    params:
        bc_threshold=10,
    shell:
        """
        mpralib combine get-reporter-variants \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} \
        --output-reporter-variants {output} > {log} 2>&1
        """


rule get_reporter_genomic_variants:
    input:
        quantification="results/{id}/quantification/{id}.{level}.variant.output.tsv.gz",
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        "results/{id}/reporter_genomic_variants/{id}.reporter_genomic_variants.{level}.bed.gz",
    log:
        "logs/variants/get_reporter_genomic_variants.{id}.{level}.log",
    benchmark:
        "benchmarks/variants/get_reporter_genomic_variants.{id}.{level}.tsv"
    wildcard_constraints:
        level="(barcode)|(oligo)",
    conda:
        getCondaEnv("mpralib.yaml")
    container:
        "docker://visze/mpralib:0.10.4"
    threads: 1
    resources:
        mem_mb=lambda wc, input: calc_mem_gb(input[1], 75) * 1024,  # Adjust memory based on input size
    params:
        bc_threshold=10,
    shell:
        """
        mpralib combine get-reporter-genomic-variants \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} \
        --output-reporter-genomic-variants >(bgzip -c > {output}) > {log} 2>&1
        """
