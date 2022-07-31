# Tool for detection of chhromothriptic rearrangements in structural variations data

---
## Methods

This tool marks with chromothripsis those breakpoints that form a spatial cluster. This assumption is based on [this](https://doi.org/10.1007/978-3-031-06220-9_13) article. In order to select from the whole set of structural variations those that form a spatial cluster, the following algorithm was used.
A weighted graph constructed, whose vertices are the breakpoints, and the weights of the edges are the values from the Hi-C matrix for the corresponding pair of the breakpoints. The resulting graph is searched for a dense subgraph, since a chromothripsis event occurs only once in one patient, one subgraph will contain all the chromothripsis breakpoints. The vertices from the resulting subgraph are translated back into breakpoints, and annotated with the structural variations in which they are included, the structural variations with both breakpoints marked as chromotripsis are labeled chromotriptic and are the final result of the algorithm.

## Installation and usage

1. Clone all repository
```bash
git clone https://github.com/Dv1t/hic_chromo_detection
```
2. Download example .mcool file [here](https://disk.yandex.ru/d/adY1-p4Nfhgj9Q) and place it in the tool folder.

3. Install required packages
```
pip install pandas, PyMaxflow, argparse, numpy, cooler
```
or
```
pip install pandas, PyMaxflow, argparse, numpy
conda install -c conda-forge -c bioconda cooler
```
4. Run main script:
```
python hic_chromo.py --sv example_data/prostata_sv.csv --hic Prostata40к80k200к400к800к.mcool --o result.csv
```

## Script parameters
Main script of project including all modules together can be run from anywhere as console tool.

### Required input

#### `--sv`
Path to the file with structural variations data, .csv extension.

#### `--hic`
Path to file with Hi-C matrix, .mcool extension. (.cool can be converted to .mcool by cooler zoomify)
[cooler-zoomify](https://cooler.readthedocs.io/en/latest/cli.html#cooler-zoomify) — how to generate a multi-resolution cooler file by coarsening.

#### `--o`
Path to output file, .csv extension.

### Optional input

#### `--r`
Hi-C matrix resolution in bases. Default is 40000.

### Output
Output file is a table in .csv format. It consist from 24 columns.
`patient_id` — ID of the patient with structural variations from the input data.
`1chr`—`22chr`,`Xchr` — arrays of coordinates of breakpoints, that algorithm marked as chromothriptic.

## Example run and data
Example structural variations data is available in `example-data` folder. Example .mcool file could be downloaded [here](https://disk.yandex.ru/d/adY1-p4Nfhgj9Q).
