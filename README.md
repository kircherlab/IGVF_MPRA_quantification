# IGVF MPRA quantification workflow

With this workflow you can use the output of MPRAsnakeflow (reporter experiment barcode file) together with the MPRA library sequence design file to quantify the MPRA reporter expression. It can quantify variants as well as elements. For elements, a label file is needed and the control and test label id of oligos. E.g. a label file looks like:

```text
shuffled_oligo_1	negative
shuffled_oligo_2	negative
test_oligo_1	test
test_oligo_2	test
test_oligo_3	test
```

Further, it uses BCalm or mpralm to quantify and make the statistics. You can control this via the config.

It generates reporter variants, reporter genomic variants, reporter elements, and reporter genomic elements file formats (defined by the IGVF MPRA FG).

## Requirements

- snakemake
- apptainer

## Variant quantification

Input files are:
- `/path/to/data/reporter_experiment_barcode.tsv.gz`: contains the counts per barcode.
- `/path/to/data/mpra_library_design.tsv.gz`: contains the MPRA library sequence design file.

Example using BCalm:

```bash
snakemake --sdm apptainer \
--config method=bcalm \
id=test \
count_file=/path/to/data/reporter_experiment_barcode.tsv.gz \
sequence_design_file=/path/to/data/mpra_library_design.tsv.gz \
--apptainer-args "-B /path/to/data -B $HOME/.cache/snakemake" \
-c 1 all_variants
```

Here we use an ID `test` so that the files are tagged with `test`. Output files are here:
 - `results/test/reporter_genomic_variants/test.reporter_genomic_variants.bcalm.tsv.gz`
 - `results/test/reporter_variants/test.reporter_variants.bcalm.tsv.gz`

With a vulcano plot:
 - `results/test/quantification/test.bcalm.variant.vulcano.png`

The config option `method=bcalm` can be removed because it is the default.

For running mpralm you can just use `method=mpralm`.

## Element quantification

Input files are:
- `/path/to/data/reporter_experiment_barcode.tsv.gz`: contains the counts per barcode.
- `/path/to/data/mpra_library_design.tsv.gz`: contains the MPRA library sequence design file.

Example using BCalm:

```bash
snakemake --sdm apptainer \
--config method=bcalm \
id=test \
count_file=/path/to/data/reporter_experiment_barcode.tsv.gz \
sequence_design_file=/path/to/data/mpra_library_design.tsv.gz \
test_label=test \
control_label=negative \
label_file=/path/to/data/labels.tsv.gz \
--apptainer-args "-B /path/to/data -B $HOME/.cache/snakemake" \
-c 1 all_elements
```

Here we use an ID `test` so that the files are tagged with `test`. Output files are here:
 - `results/test/reporter_genomic_elements/test.reporter_genomic_elements.bcalm.tsv.gz`
 - `results/test/reporter_elements/test.reporter_elements.bcalm.tsv.gz`

With a density and a vulcano plot:
 - `results/test/quantification/test.bcalm.element.density.png`
 - `results/test/quantification/test.bcalm.element.vulcano.png`

The config option `method=bcalm` can be removed because it is the default.

For running mpralm you can just use `method=mpralm`.
