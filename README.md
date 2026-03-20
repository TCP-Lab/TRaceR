# TRaceR

## Pipeline
This is a [**Kerblam!**](https://kerblam.dev/) data analysis project.
Provided that you have the right input data structure, just run the `tracer :: Multiparametric Calcium Trace Analyzer` workflow through
```bash
kerblam run tracer
```
to get the final output.

Run
```bash
kerblam run tracer 2>&1 | tee log.txt
```
to also locally save the onscreen log.

You can also run `stat_tracer :: TRaceR - Statistical Analysis Only` if you want to run only the statistical analysis step, following the time course analysis.

## Biological Model
Primary cell cultures of mouse cerebellar granules for the study of pharmacologically-induced *ataxia-telangiectasia*.

ATM (Ataxia-telangiectasia mutated) is a crucial serine/threonine protein kinase activated primarily by DNA double-strand breaks and oxidative stress.
As a central regulator of the DNA damage response (DDR), ATM phosphorylates substrates like p53, BRCA1, and Nbs1 to maintain genomic stability.
Beyond DNA damage, ATM functions as a sensor for oxidative stress, managing reactive oxygen species (ROS) levels.
Deficiencies in ATM cause Ataxia-telangiectasia, characterized by neurodegeneration.

## Experimental Conditions
- DMSO: control condition (KU vehicle - always present)
- UV: inducer of ROS production (always present)
- KU: ATM protein inhibitor (-> ataxia-telangiectasia mimicker)
- NAC: N-acetylcysteine (antioxidant agent)
- OMO: activator of NRF2 (a ROS-protective gene expression program)
- CBX: carbenoxolone (gap junction and pannexin-1 blocker)
- Gd: gadolinium (calcium channel blocker)
- TROLOX: Analogue of Vitamin E (antioxidant agent)

1. Batch **KU_lessUV_NAC**
    - DMSO
    - DMSO_NAC
    - KU
    - KU_lessUV
    - KU_NAC
1. Batch **KU_OMO**
    - DMSO_OMO
    - KU
    - KU_OMO_24h
    - KU_OMO_48h
1. Batch **KU_Gd_CBX**
    - KU
    - KU_CBX
    - KU_Gd
    - KU_Gd_CBX
1. Batch **KU_TROLOX**
    - KU
    - KU_TROLOX

## Data Format
The input data must be organized as follows: for each batch (*KU_OMO*, *KU_lessUV_NAC*, etc.), there must be a folder containing a subfolder for each condition (e.g., *DMSO*, *KU*, *KU_NAC*, etc.).
Since the names of these folders will be used to identify batches, groups, and perform comparisons, avoid using spaces in the folder names.

It is assumed that each input file—representing a biological replicate—consists of a CSV organized as follows:
 - `,` for decimal point
 - `;` as field separator
 - a string in the first row to skip
 - column headings in the second row
 - *fileEncoding* `UTF-16LE` (with BOM)

Finally, each experiment is expected to last 900 seconds, at a sampling rate of 1 second.
Otherwise, the statistics may no longer be comparable.

## Statistical Contrasts
Contrasts need to be defined inside the input folder of each batch through a plain text file, named `comp.txt` (or `comparisons.txt`), where each statistical comparison is listed in the form
```
condition -- reference
```
Blank lines allowed, as well as comment lines introduced by `#`.

## Requisites
**R packages**
- r4tcpl
- stats
- dplyr
- tidyr
- ggplot2
- patchwork
- zoo
- pracma

**Kerblam!** (>= 1.2.1)

## Anomalies
### Recording Length
The standard duration of each experiment should be exactly 900 seconds, at a sampling rate of 1 second.
The following files are exceptions.
However, given the experimental conditions, we believe they do not introduce statistical artifacts.

Only **600 time samples** for `./data/in/KU_lessUV_NAC/DMSO` experiments
- 250908_011_DMSO.csv
- 250908_012_DMSO.csv
- 250908_013_DMSO.csv

Only **700 time samples** for `./data/in/KU_lessUV_NAC/KU_NAC` experiments
- 250908_045_KU_NAC.csv

### CSV separators
`./data/in/KU_OMO/KU/251113_011_KU.csv` was the only CSV input file with a different separator convention, even though it used the same file encoding, `UTF-16LE` (with BOM).
This caused the analysis in the `./data/in/KU_OMO/KU/` folder to fail, although it did not stop the entire program from running.
It was fixed by:
```r
fpath <- "./data/in/KU_OMO/KU/251113_011_KU.csv"
fpath |> read.table(
        header = TRUE,
        sep = "\t",
        dec = ".",
        skip = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        fileEncoding = "UTF-16LE") |>
    lapply(as_num) |> as.data.frame() |> write.table(
        file = "./data/in/KU_OMO/KU/251113_011_KU_new.csv",
        sep = ";",
        dec = ",",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        fileEncoding = "UTF-16LE")
```
Then the first-line string `Channel.002\r\n` was manually reinserted in `./data/in/KU_OMO/KU/251113_011_KU_new.csv` and saved as `./data/in/KU_OMO/KU/251113_011_KU.csv`.
