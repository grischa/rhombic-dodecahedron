################################################################
#
# RHOMBIC DODECAHEDRON WATER BOX GENERATION FOR CHARMM
# TOPOLOGIES IN VMD
#
# This script generates a rhombic dodecahedron water box given
# initial CHARMM topology and coordinates. Cell dimensions and
# origins are written to stdout and saved to a file for use in NAMD.
# The input coordinates are centered in the process.
# Only for proteins!
# The larger the system, the longer the execution...
#
# To run execute: vmd -dispdev text  -e rhombicdodecahedron.tcl
#
#
# INPUT (section immediatly below this text box):
#
# psf		Input CHARMM-compatible topology.
# pdb		Input PDB coordinates.
# prefix	Prefix for segment id of water molecules.
#			See TIPS below to customize. The solvate plugin default is 'WT'.
# inc 		Parameter to speedup the script. Every 'inc' C-alpha
#			is considered in the calculation of the input molecule
#			dimension. See TIPS below to customize.
# rbddh		Output prefix for topology, coordinates and cell
#			dimensions files.
# padding   Distance between the edge of the box to the closest C-alpha.
#
#
# TIPS
#
# For small proteins (less than 1000 residues) set 'inc' to 1, it can be
# increased for very large protein systems such as entire viral
# capsids up to, say, 5.
#
# For very large proteins, the number of water segments can be such
# that a four letter limit in the segment id is breached. Set the
# variable 'prefix' to W instead of the default (WT).
#
# Created by Grischa Meyer and Cyril Reboul
#
################################################################

set psf         PSFfilename.psf
set pdb         PDBfilename.pdb
set prefix      W
set inc         5
set rbddh       Outputfilename
set padding     16.0



# DO NOT EDIT BELOW UNLESS YOU ARE AN EXPERT USER

proc dmax_alt {selection {inc 1}} {
# every 'inc' pairs calculation to work out the maximum distance between two
# atoms in an atomselection

    if {[$selection num] <= 0} {
	error "dmax: needs a selection"
    }

    set dmax2 -1.0
    set coord [$selection get {x y z}]
    set imax [expr [$selection num] - 1]
    set jmax [$selection num]
    for { set i 0 } { $i < $imax } { incr i $inc} {
	set ri [lindex $coord $i]
	for { set j [expr $i + 1]} { $j < $jmax } { incr j $inc} {
	    set rj [lindex $coord $j]
	    set rij [vecsub $ri $rj]
	    set d2 [vecdot $rij $rij]
	    if {$d2 > $dmax2} {set dmax2 $d2}
	}
    }
    return [expr sqrt($dmax2)]
}

# Main part starts here

# 'box' is only used internally for temporary files later deleted
set box         box

package require psfgen
resetpsf

mol load psf $psf pdb $pdb
set sel [atomselect top all]
set ca [atomselect top "name CA"]
puts "Calculating distance between C_alpha's. May take some time . . ."
set Dm [expr [dmax_alt $ca $inc] / 2]
set R [expr ($Dm + $padding)*sqrt(2)]
set D $R

# find mass center...
set center [lindex [measure inertia $sel] 0]
# ...and centers
$sel moveby [vecinvert $center]
set tmp tmp.pdb
$sel writepdb $tmp
mol delete all

# Solvate in big box that contains the dodecahedron and loads it
package require solvate
set xymm [expr $R*sqrt(2.)/2]
set min "-$xymm -$xymm -$R"
set max "$xymm $xymm $R"
set minmax [list $min $max]
solvate $psf $tmp -s $prefix -b 1.5 -o $box -minmax $minmax
mol delete top
resetpsf
mol load psf ${box}.psf pdb ${box}.pdb
readpsf  ${box}.psf
coordpdb ${box}.pdb

# rotate along z-axis
set all [atomselect top all]
$all move [trans z 45.0]

# Trims unnecessary waters (outside dodecahedron)
set seltext "same residue as not (x+y>-$D and x-y<$D and z-y<$D and z+y>-$D and x+z>-$D and x-z<$D and x+y<$D and z+y<$D and z+x<$D and z-y>-$D and x-z>-$D and x-y>-$D)"
set selDel [atomselect top "$seltext"]
set sel [atomselect top "not ($seltext)"]
puts "Deleting [$selDel num] atoms"
set delList [$selDel get {segid resid}]
set delList [lsort -unique $delList]

foreach record $delList {
    foreach {segid resid} $record { break }
    delatom $segid $resid
}

#rotate along z-axis again
$sel move [trans z 45.0]

$sel writepsf "$rbddh.psf"
$sel writepdb "$rbddh.pdb"

# remove temp files
file delete ${box}.log ${box}.psf ${box}.pdb combine.pdb combine.psf $tmp temp.psf temp.pdb

# Some info
set D [expr $D*sqrt(2)]
puts "Center of mass is at: 0.0 0.0 0.0"
puts "Radius: $D"
puts "Maximum distance between C-alphas: [expr $Dm * 2.0]"
puts "Dimensions:"
set a "cellBasisVector1 $D 0.0 0.0"
set b "cellBasisVector2 0.0 $D 0.0"
set c "cellBasisVector3 [expr $D/2] [expr $D/2] [expr $D/sqrt(2)]"
set fh [open "${rbddh}_cell_basis_vectors.text" "w"]
puts $fh "cellBasisVector1 $D 0.0 0.0"
puts $fh "cellBasisVector2 0.0 $D 0.0"
puts $fh  "cellBasisVector3 [expr $D/2] [expr $D/2] [expr $D/sqrt(2)]"
close $fh
exit
