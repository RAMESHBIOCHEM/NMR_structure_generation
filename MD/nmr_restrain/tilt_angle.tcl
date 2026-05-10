
#close output
#close $outdiid
mol new ../step5_input.psf type psf
mol addfile ../step5_input.pdb
#mol addfile ../step5_input.pdb type pdb first 0 last -1 step 1 waitfor all
set all [atomselect top "all"]
set PROA [atomselect top "protein and segname PROA and name CA"]
set PROB [atomselect top "protein and segname PROB and name CA"]
$all set occupancy 1.0
$all set beta 0.0
$PROA set beta 1.0
$all writepdb PROA.ref

$all set beta 0.0
$PROB set beta 1.0
$all writepdb PROB.ref

$all delete
$PROA delete
$PROB delete

#}
quit
