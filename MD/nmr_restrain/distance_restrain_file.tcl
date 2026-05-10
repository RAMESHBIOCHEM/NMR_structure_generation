
#close output
#close $outdiid
mol new ../step5_input.psf type psf
mol addfile ../step5_input.pdb
#mol addfile ../step5_input.pdb type pdb first 0 last -1 step 1 waitfor all
set outfile [open Distance.ref w]
##3881

#set j [expr {$i+1}]
#set k [expr $i+2}]

#set k [expr {$i+2}]
set i1 [atomselect top "segname PROA and name CA and resid 17"]
set j1 [atomselect top "segname PROB and name CA and resid 18"]
set Ai1 [$i1 get index]
set Bj1 [$j1 get index]
set bond_AB1 [list $Ai1 $Bj1]

set i2 [atomselect top "segname PROA and name CA and resid 28"]
set j2 [atomselect top "segname PROB and name CA and resid 29"]
set Ai2 [$i2 get index]
set Bj2 [$j2 get index]
set bond_AB2 [list $Ai2 $Bj2]

set i3 [atomselect top "segname PROA and name CA and resid 25"]
set j3 [atomselect top "segname PROB and name CA and resid 25"]
set Ai3 [$i3 get index]
set Bj3 [$j3 get index]
set bond_AB3 [list $Ai3 $Bj3]


   # set rlist5 [$PS2 get index]
    puts $outfile "bond $bond_AB1 ZFC 6.0\nbond $bond_AB2 ZFC 6.0\nbond $bond_AB3 ZFC 5.38"
## grep -v "DIHEDRAL  ZFC 0.0" dihe.txt > temp.txt
## grep -v "IMPROPER   ZFC 120.0" temp.txt > dihe.txt
## open vi editor for dihe.txt, then type :%s/ZFC/$FC/g
    close $outfile
#}
quit
