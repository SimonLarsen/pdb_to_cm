pdb\_to\_cm
===========

Simple Python script for computing protein contact maps from PDB files as described in Godzik and Skolnick, (1994).
An edge is added between all non-adjacent pairs of amino acids where the euclidean distance between their alpha carbons are less then some threshold.

## Usage

Compute the contact map for 1BPI with a threshold of 7.5 ångström.

```bash
python pdb_to_cm.py 1bpi.pdb 1bpi.cm -t 7.5
```

## References

* Godzik A, Skolnick J. Flexible algorithm for direct multiple alignment of protein structures and sequences. Computer applications in the biosciences: CABIOS. 1994 Dec 1;10(6):587-96.
