# prot-fin experimental space

You find the reference proteins [here](https://github.com/usadellab/prot-fin/raw/5be77c4247327e3958c89200c03a938ec4734834/material/Mapman_reference_DB_202310.tar.bz2)

The table of Kidera factors is from [here](https://github.com/usadellab/prot-fin/raw/5be77c4247327e3958c89200c03a938ec4734834/material/Amino_Acid_Kidera_Factors.csv), the result of 

## prot-fin

Usadellab's contributions to the prot-fin project.

### Project location on the IBG-4 cluster

Go to `/mnt/data/hakimeh/prot-fin` to find the local clone of this repository.

### Project structure

The following directories contain the content indicated by their names:
- `docs`
- `experiments`
- `materials`

### Generate table of Kidera Factors

Kidera Factors are numeric values that describe the physical and chemical
properties of amino acids, e.g. hydrophobicity or volume. Oversimplified they
are derived from a principal component analysis of more than 180 physical and
chemical features of amino acids. For the original reference see below or the R
package `Peptides`.

Use [this R-script](https://github.com/usadellab/prot-fin/blob/5be77c4247327e3958c89200c03a938ec4734834/methods/Amino_Acid_Kidera_Factors.R) to generate output table
`./materials/Amino_Acid_Kidera_Factors.csv`.

On Kidera Factors please see:
> Kidera, A., Konishi, Y., Oka, M., Ooi, T., & Scheraga, H. A. (1985).
> Statistical analysis of the physical properties of the 20 naturally occurring
> amino acids. Journal of Protein Chemistry, 4(1), 23-55.

### Reference proteins

We use the reference proteins (amino acid sequences in Fasta format) that were
used to generate the Mapman Bin hidden Markov Models (HMMs). The two files,
i.e. the Fasta file with the amino acid sequences and their UniProt identifiers
and the file annotating these reference proteins with the Mapman Bin they
belong to are compressed in the archive
[`Mapman_reference_DB_202310.tar.bz2`](https://github.com/usadellab/prot-fin/raw/5be77c4247327e3958c89200c03a938ec4734834/material/Mapman_reference_DB_202310.tar.bz2). It contains the two files:

- `mapmanreferencebins.results.txt` - The file assigning Mapman Bins to
  reference proteins.
- `protein.fa` - The Fasta file with the UniProt identifiers and their amino
  acid sequences.

Note that the uncompressed content of the above archive is ignored by `git` to avoid big data issues.