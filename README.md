Rhombic Dodecahedron Water Box Generation For Charmm Topologies In VMD
======================================================================

This script generates a rhombic dodecahedron water box given
initial CHARMM topology and coordinates. Cell dimensions and
origins are written to stdout and saved to a file for use in NAMD.

The input coordinates are centered in the process.

Only for proteins!

The larger the system, the longer the execution...

To run execute: `vmd -dispdev text -e rhombicdodecahedron.tcl`

## Input

psf    | Input CHARMM-compatible topology.
pdb    | Input PDB coordinates.
prefix | Prefix for segment id of water molecules.
       | See TIPS below to customize. The solvate plugin default is 'WT'.
inc    | Parameter to speedup the script. Every 'inc' C-alpha
       | is considered in the calculation of the input molecule
       | dimension. See TIPS below to customize.
rbddh  | Output prefix for topology, coordinates and cell
       | dimensions files.
padding| Distance between the edge of the box to the closest C-alpha.

## Tips

For small proteins (less than 1000 residues) set 'inc' to 1, it can be
increased for very large protein systems such as entire viral
capsids up to, say, 5.

For very large proteins, the number of water segments can be such
that a four letter limit in the segment id is breached. Set the
variable 'prefix' to W instead of the default (WT).

Created by Grischa Meyer and Cyril Reboul
