
#close output
#close $outdiid
mol new ../step5_input.psf type psf
mol addfile ../step5_input.pdb
#mol addfile ../step5_input.pdb type pdb first 0 last -1 step 1 waitfor all
set outfile [open Dihedral3.ref w]
##3881

#set j [expr {$i+1}]
#set k [expr $i+2}]

for {set i 12} {$i < 36} {incr i} {

set j [expr {$i+1}]
set k [expr {$i+2}]
set i_C [atomselect top "segname PROA and name C and resid $i"]
set j_N [atomselect top "segname PROA and name N and resid $j"]
set j_CA [atomselect top "segname PROA and name CA and resid $j"]
set j_C [atomselect top "segname PROA and name C and resid $j"]

set Ci [$i_C get index]
set Nj [$j_N get index]
set CAj [$j_CA get index]
set Cj [$j_C get index]
set phi [list $Ci $Nj $CAj $Cj]

####################################################################
set k_N [atomselect top "segname PROA and name N and resid $k"]
set Nk [$k_N get index]
set psi [list $Nj $CAj $Cj $Nk]

   # set rlist4 [$PS1 get index]
   # set rlist5 [$PS2 get index]
    puts $outfile "dihedral $phi ZFC -60.0 \ndihedral $psi ZFC -45.0"
    }
 
   # puts $outfile "*************************************"
for {set i 12} {$i < 36} {incr i} {

set j [expr {$i+1}]
set k [expr {$i+2}]
set i_C [atomselect top "segname PROB and name C and resid $i"]
set j_N [atomselect top "segname PROB and name N and resid $j"]
set j_CA [atomselect top "segname PROB and name CA and resid $j"]
set j_C [atomselect top "segname PROB and name C and resid $j"]

set Ci [$i_C get index]
set Nj [$j_N get index]
set CAj [$j_CA get index]
set Cj [$j_C get index]
set phi [list $Ci $Nj $CAj $Cj]

####################################################################
set k_N [atomselect top "segname PROB and name N and resid $k"]
set Nk [$k_N get index]
set psi [list $Nj $CAj $Cj $Nk]

   # set rlist4 [$PS1 get index]
   # set rlist5 [$PS2 get index]
    puts $outfile "dihedral $phi ZFC -60.0 \ndihedral $psi ZFC -45.0"
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
