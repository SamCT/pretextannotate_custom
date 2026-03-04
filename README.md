# pretextannotate

PretextAnnotate is a script originally written by Karen van Niekerk (Sanger GRIT) in order to add chromosome annotations to a PNG image produced by PretextSnapshot.

`/src/fonts/OpenSans-Regular.ttf` - Open Sans font file used for text rendering.
Font file is taken from https://fonts.google.com/specimen/Open+Sans/ and is licensed under [Open Font License](https://openfontlicense.org/).

`original_scripts` - Contains the original scripts used in this project, this is simply for archiving and reference purposes.

## Installation
```
git clone https://github.com/sanger-tol/pretextannotate.git

cd pretextannotate

pip install .
```

#### Coming soon
```
pip install pretextannotate
```

## Usage

To make full use of the script (and politeness to NCBI) you must set the environment variables:
`ENTREZ_EMAIL`
`ENTREZ_API_KEY`

You can generate your own API keys by creating an account on `https://www.ncbi.nlm.nih.gov/account/` with your ORCID.

```
pretextannotate -h

pretextannotate \\
    --pretext_file src/tests/ilDryDodo1.1_normal_FullMap.png \\
    --output ./ \\
    --prefix HELLO \\
    --gca_accession GCA_965178025.1
    
    
pretextannotate \\
    --pretext_file src/tests/ilDryDodo1.1_normal_FullMap.png \\
    --output ./ \\
    --prefix HELLO \\
    --sizes GCA_965178025.1.sizes
```

### Custom context files (no GCA accession required)

If you do not have a `GCA_*` accession yet, you can build a context file directly from your FASTA headers (for example `>H1.chr01`, `>H2.chr01`, ...):

```bash
pretextannotate-build-context \
    --fasta my_assembly.fa \
    --mapping molecule_mapping.tsv \
    --output custom_context.tsv
```

Then run `pretextannotate` with that file:

```bash
pretextannotate \
    --pretext_file my_pretext.png \
    --sizes custom_context.tsv \
    --vertical_label_field INSDC \
    --prefix MY_SAMPLE \
    --output ./
```

Example `molecule_mapping.tsv` (optional):

```tsv
# sequence_name	molecule_label
H1.chr01	1
H1.chr02	2
H2.chr01	3
H2.chr02	4
```

Example generated `custom_context.tsv`:

```tsv
# INSDC	length_bp	molecule
H1.chr01	50324412	1
H1.chr02	48219877	2
H2.chr01	49910332	3
H2.chr02	48000741	4
```

## Expected output

With the input pretext snapshot:
![ilDryDodo1 pretext map](./src/tests/ilDryDodo1.1_normal_FullMap.png)

As well as the arguments used in the Usage section.

The output should be a PNG, gif and tif file resembling:
![ilDryDodo1 annotatedpretext map](./src/tests/ilDryDodo1_Pretext.png)

There will also be a pretextannotation.log file containing, in this case:
```
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Starting Pretext Annotation
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] PretextSnapshot: src/tests/ilDryDodo1.1_normal_FullMap.png | WITH | context_dict: {"accession": "GCA_965178025.1"}
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Input Snapshot Image is src/tests/ilDryDodo1.1_normal_FullMap.png
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Output will be saved at .//HELLO_annotated_pretext.png
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Starting Pretext Annotation Process
2026-02-12 12:55:46,470 [INFO] [Pretext Annotation] Adjusted font size: 60 for 31 chromosomes
2026-02-12 12:55:47,355 [INFO] [Pretext Annotation] Saved labelled PNG → .//HELLO_annotated_pretext.png
2026-02-12 12:55:47,783 [INFO] [Pretext Annotation] Converted .//HELLO_annotated_pretext.png → .//HELLO_annotated_pretext.tif, .//HELLO_annotated_pretext.gif
2026-02-12 12:55:47,783 [INFO] [Pretext Annotation] Converted to TIFF & GIF → .//HELLO_annotated_pretext.tif, .//HELLO_annotated_pretext.gif
```

## Future ToDo's:
- Graphs to right side
    - Telomere, gap, coverage, repeats, GC?
- Tests
