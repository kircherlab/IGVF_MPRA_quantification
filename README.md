# IGVF MPRA quantification workflow

With this workflow you can use the output of MPRAsnakeflow (reporter experiment barcode file) together with the MPRA library sequence design file to quantify MPRA reporter expression. It can quantify variants as well as elements. For elements, a label file is needed and the control and test label IDs of oligos. For example, a label file looks like:

```text
shuffled_oligo_1	negative
shuffled_oligo_2	negative
test_oligo_1	test
test_oligo_2	test
test_oligo_3	test
```

Further, it can quantify using all available barcodes per oligo (config `level=barcode`) or it can use an aggregated version (config `level=oligo`). For quantification and statistics it uses BCalm.

```text
Keukeleire, Pia, Jonathan D. Rosen, Angelina Göbel-Knapp, Kilian Salomon, Max Schubach, and Martin Kircher. 
"Using individual barcodes to increase quantification power of massively parallel reporter assays." 
BMC Bioinformatics 26, no. 1 (2025): 52.
```

It generates reporter variants, reporter genomic variants, reporter elements, and reporter genomic elements file formats (defined by the IGVF MPRA FG).

## Requirements

- snakemake
- Apptainer/Singularity

## Variant quantification

Input files are:
- `/path/to/data/reporter_experiment_barcode.tsv.gz`: contains the counts per barcode.
- `/path/to/data/mpra_library_design.tsv.gz`: contains the MPRA library sequence design file.

Example using barcode-level quantification:

```bash
snakemake --sdm apptainer \
--config level=barcode \
id=test \
count_file=/path/to/data/reporter_experiment_barcode.tsv.gz \
sequence_design_file=/path/to/data/mpra_library_design.tsv.gz \
--apptainer-args "-B /path/to/data -B $HOME/.cache/snakemake" \
-c 1 all_variants
```

Here we use an ID `test` so that the files are tagged with `test`. Output files are here:
 - `results/test/reporter_genomic_variants/test.reporter_genomic_variants.barcode.tsv.gz`
 - `results/test/reporter_variants/test.reporter_variants.barcode.tsv.gz`

With a volcano plot:
 - `results/test/quantification/test.barcode.variant.volcano.png`

The config option `level=barcode` can be removed because it is the default.

For running on aggregated counts you can just use `level=oligo`.

## Element quantification

Input files are:
- `/path/to/data/reporter_experiment_barcode.tsv.gz`: contains the counts per barcode.
- `/path/to/data/mpra_library_design.tsv.gz`: contains the MPRA library sequence design file.

Example using barcode-level quantification:

```bash
snakemake --sdm apptainer \
--config level=barcode \
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
 - `results/test/reporter_genomic_elements/test.reporter_genomic_elements.barcode.tsv.gz`
 - `results/test/reporter_elements/test.reporter_elements.barcode.tsv.gz`

With a density and a volcano plot:
 - `results/test/quantification/test.barcode.element.density.png`
 - `results/test/quantification/test.barcode.element.volcano.png`

The config option `level=barcode` can be removed because it is the default.

For running on aggregated counts you can just use `level=oligo`.
