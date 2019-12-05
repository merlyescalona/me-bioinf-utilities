# BioNano discordant regions - HiRise filtering

## Installation 

Tested on:
- R 3.4.4 "Someone to lean on".
- Platform: x86_64-pc-linux-gnu (64-bit)

### Requirements 
- rlang (0.2.1)
- stringi (1.4.3)
- stringr (1.3.1)
- gdata (2.18.0)
- egg (0.4.5)
- optparse (1.6.4)

## Usage

```
Rscript --vanilla discordant_regions_filter.R -d regions.bed -g gaps.bed -s scaffold.sizes [-t THRESHOLD -p PREFIX -o OUTPUT_FOLDER]
```

## Options:
- `-d CHARACTER, --discordant=CHARACTER`: Path of the BED file describing the discordant regions.

- `-s CHARACTER, --sizes=CHARACTER`: Path of the TSV file with the scaffold sizes.
    - Format: 
    ```
    scaffold    size
    ```
- `-g CHARACTER, --gaps=CHARACTER`: Path of the BED file describing the gaps in the scaffolds.
    - Format: 
    ```
    scaffold   start   end
    ```
- `-t NUMBER, --threshold=NUMBER`: Distance threshold to define closeness of a discordant region to a gap [default= 50,000]

- `-p CHARACTER, --prefix=CHARACTER`: Output prefix  [default= out]
- `-o CHARACTER, --output=CHARACTER`: Path to the output folder [default= ~/]
- `-h, --help`: Show this help message and exit

## Output
- This Rscript generates a list of BED files:
    - `prefix.breakpoints.bed`: 
    Breakpoints file to be used directly to update the genome.
        - Format: `scaffold   start   end     -`
        - Where "-" defines that the specified gap/join must be broken.
    - `prefix.hirise.auto.bed: `: 
    - `prefix.hirise.auto.conflicts.bed`:
    - `prefix.hirise.auto.conflicts_inv_invpart.bed`:
    - `prefix.hirise.auto.conflicts_invpair.bed`:
    - `prefix.hirise.auto.conflicts_invpart.bed`:
    - `prefix.hirise.auto.inv.bed`:
    - `prefix.hirise.auto.inv_invpart.bed`:
    - `prefix.hirise.auto.invpair.bed`:
    - `prefix.hirise.auto.invpair_invpart.bed`:
    - `prefix.hirise.auto.invpart.bed`:

## Wiki

## Troubleshooting/bugs/question

- Please create an issue.

## Contribution 

(c) 2019 Merly Escalona <mmescalo@ucsc.edu>

Paleogenomics Lab.

Department of Biomolecular Engieering.

University of California Santa Cruz

