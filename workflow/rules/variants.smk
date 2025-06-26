rule get_variant_counts:
    # container:
    #   "docker://quay.io/biocontainers/mpralib:0.7.3--pyhdfd78af_0"
    conda:
        "mpralib"
    input:
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        variant_counts="results/{id}/quantification/{id}.{method}.variant.input.tsv.gz",
    log:
        "logs/variants/get_variant_counts.{id}.{method}.log",
    params:
        normalize="--normalized-counts" if config["mpralib_normalized_counts"] else "",
        bc_threshold=1,
        barcodes=lambda wc: "--barcodes" if wc.method == "bcalm" else "--oligos",
    shell:
        """
        mpralib sequence-design get-variant-counts \
        --input {input.counts} --sequence-design {input.sequence_design} \
        {params.barcodes} --bc-threshold {params.bc_threshold} {params.normalize} \
        --output {output.variant_counts} > {log} 2>&1
        """


rule get_variant_map:
    # container:
    #     "docker://quay.io/biocontainers/mpralib:0.7.3--pyhdfd78af_0"
    conda:
        "mpralib"
    input:
        sequence_design=config["sequence_design_file"],
    output:
        variant_map="results/{id}/{id}.variant_map.tsv.gz",
    log:
        "logs/variants/get_variaget_variant_mapnt_counts.{id}.log",
    shell:
        """
        mpralib sequence-design get-variant-map \
        --sequence-design {input.sequence_design} \
        --output {output.variant_map} > {log} 2>&1
        """


rule run_variants_bcalm_quantification:
    container:
        "docker://visze/bcalm:latest"
    input:
        variant_counts="results/{id}/quantification/{id}.bcalm.variant.input.tsv.gz",
        variant_map="results/{id}/{id}.variant_map.tsv.gz",
        script=getScript("bcalm_variants.R"),
    output:
        result="results/{id}/quantification/{id}.bcalm.variant.output.tsv.gz",
        vulcano_plot="results/{id}/quantification/{id}.bcalm.variant.vulcano.png",
    log:
        "logs/variants/run_variants_bcalm_quantification.{id}.log",
    params:
        normalize="FALSE" if config["mpralib_normalized_counts"] else "TRUE",
    shell:
        """
        Rscript {input.script} \
        --count {input.variant_counts} --map {input.variant_map} \
        --normalize {params.normalize} \
        --output {output.result} --output-plot {output.vulcano_plot} > {log} 2>&1
        """

rule run_variants_mpralm_quantification:
    container:
        "docker://visze/bcalm:latest"
    input:
        variant_counts="results/{id}/quantification/{id}.mpralm.variant.input.tsv.gz",
        script=getScript("mpralm_variants.R"),
    output:
        result="results/{id}/quantification/{id}.mpralm.variant.output.tsv.gz",
        vulcano_plot="results/{id}/quantification/{id}.mpralm.variant.vulcano.png",
    log:
        "logs/variants/run_variants_mpralm_quantification.{id}.log",
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
    # container:
    #     "docker://quay.io/biocontainers/mpralib:0.7.3--pyhdfd78af_0"
    conda:
        "mpralib"
    input:
        quantification="results/{id}/quantification/{id}.{method}.variant.output.tsv.gz",
        counts=config["count_file"],
        sequence_design=config["sequence_design_file"],
    output:
        reporter_variants="results/{id}/{format}/{id}.{format}.{method}.tsv.gz",
    log:
        "logs/variants/get_reporter_variants.{id}.{method}.{format}log",
    wildcard_constraints:
        format="(reporter_variants)|(reporter_genomic_variants)",
    params:
        bc_threshold=10,
        output_format=lambda wc: (
            ("get-reporter-variants", "--output-reporter-variants")
            if wc.format == "reporter_variants"
            else (
                "get-reporter-genomic-variants",
                "--output-reporter-genomic-variants",
            )
        ),
    shell:
        """
        mpralib sequence-design {params.output_format[0]} \
        --input {input.counts} \
        --sequence-design {input.sequence_design} \
        --bc-threshold {params.bc_threshold} \
        --statistics {input.quantification} \
        {params.output_format[1]} {output.reporter_variants} > {log} 2>&1
        """
