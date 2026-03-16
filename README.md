# TRaceR

ATM (Ataxia-telangiectasia mutated) is a crucial serine/threonine protein kinase activated primarily by DNA double-strand breaks and oxidative stress. As a central regulator of the DNA damage response (DDR), ATM phosphorylates substrates like p53, BRCA1, and Nbs1 to maintain genomic stability. Beyond DNA damage, ATM functions as a sensor for oxidative stress, managing reactive oxygen species (ROS) levels. Deficiencies in ATM cause Ataxia-telangiectasia, characterized by neurodegeneration.

mouse cerebellar granules

DMSO: control condition (KU vehicle) (always present)
UV: inducer of ROS production (always present)
KU: ATM protein inhibitor (-> ataxia-telangiectasia)
NAC: N-acetylcysteine (antioxidant agent)
OMO: NRF2 activator (ROS-protective gene expression program)
CBX: carbenoxolone (?)
Gd: gadolinium (?)
TROLOX: (?)

Folders KU lessUV NAC
 - DMSO
 - DMSO_NAC
 - KU
 - KU_lessUV
 - KU_NAC

Folders KU OMO
 - DMSO_OMO
 - KU
 - KU_OMO_24h
 - KU_OMO_48h

Folders KU Gd CBX
 - KU
 - KU_CBX
 - KU_Gd
 - KU_Gd_CBX

Folders KU TROLOX
 - KU
 - KU_TROLOX


Comparisons are define in comps.txt files
-- add a check for groups

## Requisites
The input data must be organized as follows: for each batch (*KU_OMO*, *KU_lessUV_NAC*, etc.), there must be a folder containing a subfolder for each condition (e.g., *DMSO*, *KU*, *KU_NAC*, etc.).
Since the names of these folders will be used to identify batches, groups and perform comparisons, avoid using spaces in the folder names.

It is assumed that each input file, representing a biological replicate, consists of a CSV organized as follows:
 - `,` for decimal point
 - `;` as field separator
 - a string in the first row to skip
 - column headings in the second row
 - *fileEncoding* `UTF-16LE` (with BOM)

Finally, each experiment is expected to last 900 seconds at a sampling rate of 1 second.
Otherwise, the statistics will no longer be comparable.



## Anomalies
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
