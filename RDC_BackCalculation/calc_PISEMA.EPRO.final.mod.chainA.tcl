source /home/ramesh/nmda_receptor/nmda.ramesh/tools/bigdcd.tcl
#package require bigdcd

# rotate vector v around u
proc vecrot {u v} {
  global sin_delta cos_delta

  foreach {ux uy uz} $u {break}
  foreach {vx vy vz} $v {break}

  set oneMinusCosDelta [expr {1. - $cos_delta}]
  set R11 [expr {$cos_delta + $ux * $ux * $oneMinusCosDelta}]
  set R12 [expr {$ux * $uy * $oneMinusCosDelta - $uz * $sin_delta}]
  set R13 [expr {$ux * $uz * $oneMinusCosDelta + $uy * $sin_delta}]
  set R21 [expr {$uy * $ux * $oneMinusCosDelta + $uz * $sin_delta}]
  set R22 [expr {$cos_delta + $uy * $uy * $oneMinusCosDelta}]
  set R23 [expr {$uy * $uz * $oneMinusCosDelta - $ux * $sin_delta}]
  set R31 [expr {$uz * $ux * $oneMinusCosDelta - $uy * $sin_delta}]
  set R32 [expr {$uz * $uy * $oneMinusCosDelta + $ux * $sin_delta}]
  set R33 [expr {$cos_delta + $uz * $uz * $oneMinusCosDelta}]

  set Rvx [expr {$R11 * $vx + $R12 * $vy + $R13 * $vz}]
  set Rvy [expr {$R21 * $vx + $R22 * $vy + $R23 * $vz}]
  set Rvz [expr {$R31 * $vx + $R32 * $vy + $R33 * $vz}]

  set Rv [list $Rvx $Rvy $Rvz]
  return $Rv
}

proc calc_PISEMA {frame} {
  global num_frames
  global DC_sim_list CS_sim_list resids_CS
  global nu0 sigma11 sigma22 sigma33
    
  foreach res $resids_CS {
    set atomC [atomselect top "segname PROA and resid [expr {$res-1}] and name C"]
    set atomN [atomselect top "segname PROA and resid $res and name N"]
    set atomH [atomselect top "segname PROA and resid $res and name HN"]
    set vecC [list [$atomC get x] [$atomC get y] [$atomC get z]]
    set vecN [list [$atomN get x] [$atomN get y] [$atomN get z]]
    set vecH [list [$atomH get x] [$atomH get y] [$atomH get z]]
    set vecNH [vecnorm [vecsub $vecH $vecN]]
    set vecNC [vecnorm [vecsub $vecC $vecN]]
    
    # dipolar coupling
    set zNH [lindex $vecNH 2]
    set DC_sim [expr {0.5 * $nu0 * (3.0 * $zNH * $zNH - 1.0)}]
    set DC_sim_list($res) [expr {$DC_sim_list($res) + $DC_sim}]
    
    # chemical shift
    set e2 [vecnorm [veccross $vecNC $vecNH]]
    set e3 [vecnorm [vecrot $e2 $vecNH]]
    set e1 [vecnorm [veccross $e2 $e3]]
    set e1z [lindex $e1 2]
    set e2z [lindex $e2 2]
    set e3z [lindex $e3 2]

    set CS_sim [expr {$sigma11 * $e1z * $e1z + $sigma22 * $e2z * $e2z + $sigma33 * $e3z * $e3z}]
    set CS_sim_list($res) [expr {$CS_sim_list($res) + $CS_sim}]
    
    $atomC delete
    $atomN delete
    $atomH delete
  }
  incr num_frames
}

set num_frames 0
set exp_DC_file dc_exp.xvg
set exp_CS_file cs_exp.xvg

set sigma11 57.3
set sigma22 81.2
set sigma33 227.8
#set sigma33A 232.1
set nu0 10.735
set DEG2RAD 0.01745329251994329577
set delta -17.
set sin_delta [expr {sin($delta * $DEG2RAD)}]
set cos_delta [expr {cos($delta * $DEG2RAD)}]

set resids_CS []

#set resids_CS 13
#lappend resids_CS 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37


set inStream [open $exp_CS_file r]
while {[gets $inStream line] > 0} {
  lappend resids_CS [lindex $line 0]
}
close $inStream

array set DC_sim_list {}
foreach res $resids_CS {
  set DC_sim_list($res) 0.
}

array set CS_sim_list {}
foreach res $resids_CS {
  set CS_sim_list($res) 0.
}

mol load psf /home/ramesh/XPLOR-NIH/PDB_structure/E_peptide/traj/Epro_dimer/namd_pisema_restraint/step5_input.psf
bigdcd calc_PISEMA /home/ramesh/XPLOR-NIH/PDB_structure/E_peptide/traj/Epro_dimer/namd_pisema_restraint/bcset1/step8.1_prod.dcd
bigdcd_wait

set outStream [open dc_sim.chainA.xvg w]
foreach res $resids_CS {
  set DC_sim_list($res) [expr {$DC_sim_list($res) / double($num_frames)}]
  puts $outStream [format "%i %8.3f" $res $DC_sim_list($res)]
}
close $outStream

set outStream [open cs_sim.chainA.xvg w]
foreach res $resids_CS {
  set CS_sim_list($res) [expr {$CS_sim_list($res) / double($num_frames)}]
  puts $outStream [format "%i %8.3f" $res $CS_sim_list($res)]
}
close $outStream

exit




#mol load psf ../../../03prepProtLip/prot_lip_sol_ions.psf
#bigdcd calc_tilt ../../traj/md.1.dcd




