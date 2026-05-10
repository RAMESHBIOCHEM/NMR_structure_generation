
#close output
#close $outdiid
mol new ../step5_input.psf type psf
mol addfile ../step5_input.pdb
#mol addfile ../step5_input.pdb type pdb first 0 last -1 step 1 waitfor all
set outfile [open Hbond3.ref w]
##3881

#set j [expr {$i+1}]
#set k [expr $i+2}]

for {set i 12} {$i < 34} {incr i} {

set j [expr {$i+4}]
#set k [expr {$i+2}]
set i_O [atomselect top "segname PROA and name O and resid $i"]
set j_HN [atomselect top "segname PROA and name HN and resid $j"]


set Oi [$i_O get index]
set HNj [$j_HN get index]
set bond_OH [list $Oi $HNj]

set j_N [atomselect top "segname PROA and name N and resid $j"]
set Nj [$j_N get index]
set bond_ON [list $Oi $Nj]

   # set rlist5 [$PS2 get index]
    puts $outfile "bond $bond_OH ZFC 2.0\nbond $bond_ON ZFC 3.0"
    }
 
   # puts $outfile "*************************************"
for {set i 12} {$i < 34} {incr i} {

set j [expr {$i+4}]
#set k [expr {$i+2}]
set i_O [atomselect top "segname PROB and name O and resid $i"]
set j_HN [atomselect top "segname PROB and name HN and resid $j"]

set Oi [$i_O get index]
set HNj [$j_HN get index]
set bond_OH [list $Oi $HNj]

set j_N [atomselect top "segname PROB and name N and resid $j"]
set Nj [$j_N get index]
set bond_ON [list $Oi $Nj]
   # set rlist5 [$PS2 get index]
    puts $outfile "bond $bond_OH ZFC 2.0\nbond $bond_ON ZFC 3.0"

    }
## grep -v "DIHEDRAL  ZFC 0.0" dihe.txt > temp.txt
## grep -v "IMPROPER   ZFC 120.0" temp.txt > dihe.txt
## open vi editor for dihe.txt, then type :%s/ZFC/$FC/g
    close $outfile



#puts $outfile " DIHEDRAL $rlist1 ZFC 0.0 \nDIHEDRAL $rlist2 ZFC 0.0 \nIMPROPER $rlist3 $rlist4 ZFC 120.0"


#proc center {molid} {
#set sel1 [atomselect top "protein and resid 354"]
#set nf [molinfo top get numframes]
#set outfile [open rr.dat w]
#for {set j 1} {$j < 10} {incr j} {
#set sel2 [atomselect top "resname BGLC and resid $j"]
#  for {set i 0} {$i<$nf} {incr i} {
#puts "frame $i of $nf"
#$sel1 frame $i
#$sel2 frame $i
#set com1 [measure center $sel1 weight mass]
#set com2 [measure center $sel2 weight mass]
#set xc [lindex $com2 0]
#set yc [lindex $com2 1]
#set zc [lindex $com2 2]
#set sel3 [atomselect top "resname BGLC and ((x-$xc)*(x-$xc) + (y-$yc)*(y-$yc) + (z-$zc)*(z-$zc)) < 1.3*1.3" frame $i]
#        set rlist [$sel3 get name]
#        set slist [$sel3 get serial]
       

  #     puts $outfile "$i $j $rlist $slist"
#}
#}
#close $outfile
#}
quit
