set sigma11 57.3
set sigma22 81.2
set sigma33 227.8
set nu0 10.735
#  angle between NH vector and e3 vector
set delta -17.
set DEG2RAD 0.01745
set sin_delta [expr {sin($delta*$DEG2RAD)}]
set cos_delta [expr {cos($delta*$DEG2RAD)}]

# eq. s4
set d_xNH_d_xN -1.
set d_xNH_d_yN 0.
set d_xNH_d_zN 0.
set d_yNH_d_xN 0.
set d_yNH_d_yN -1.
set d_yNH_d_zN 0.
set d_zNH_d_xN 0.
set d_zNH_d_yN 0.
set d_zNH_d_zN -1.
# eq. s5
set d_xNH_d_xH 1.
set d_xNH_d_yH 0.
set d_xNH_d_zH 0.
set d_yNH_d_xH 0.
set d_yNH_d_yH 1.
set d_yNH_d_zH 0.
set d_zNH_d_xH 0.
set d_zNH_d_yH 0.
set d_zNH_d_zH 1.
# eq. s13
set d_xNC_d_xN -1.
set d_xNC_d_yN 0.
set d_xNC_d_zN 0.
set d_yNC_d_xN 0.
set d_yNC_d_yN -1.
set d_yNC_d_zN 0.
set d_zNC_d_xN 0.
set d_zNC_d_yN 0.
set d_zNC_d_zN -1.
# eq. s14
set d_xNC_d_xC 1.
set d_xNC_d_yC 0.
set d_xNC_d_zC 0.
set d_yNC_d_xC 0.
set d_yNC_d_yC 1.
set d_yNC_d_zC 0.
set d_zNC_d_xC 0.
set d_zNC_d_yC 0.
set d_zNC_d_zC 1.

# implict derivatives
set d_xNH_d_xC 0.
set d_xNH_d_yC 0.
set d_xNH_d_zC 0.
set d_yNH_d_xC 0.
set d_yNH_d_yC 0.
set d_yNH_d_zC 0.
set d_zNH_d_xC 0.
set d_zNH_d_yC 0.
set d_zNH_d_zC 0.

proc vecnorm {vec} {
  return [vecscale [expr {1./[veclength $vec]}] $vec]
}

proc max { a b} {
	if {$a<$b} {
		return $b
		} else {
			return $a}
	}

# according to wiki, google docs: solid state nmr
proc veccross {u v} {
  foreach {u1 u2 u3} $u {break}
  foreach {v1 v2 v3} $v {break}
  set cross1 [expr {$u2*$v3 - $u3*$v2}]
  set cross2 [expr {$u3*$v1 - $u1*$v3}]
  set cross3 [expr {$u1*$v2 - $u2*$v1}]
  return [list $cross1 $cross2 $cross3]
}

# R is from wiki. u is the normalized axis
# v is the vector that needs to be rotated
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

set DCs_exp []
set CSs_exp []
# only the residue IDs of CS need to be recorded because
# DC will also use those atomIDs
set resids_CS []
set atomIDs_CS []

set inStream_DC [open $exp_DC_file r]
while {[gets $inStream_DC line] > 0} {
  lappend DCs_exp [lindex $line 1]
}
close $inStream_DC

set inStream_CS [open $exp_CS_file r]
while {[gets $inStream_CS line] > 0} {
  lappend resids_CS [lindex $line 0]
#  puts "$resids_CS"
  lappend CSs_exp [lindex $line 1]
}
close $inStream_CS

## put last residue number of PROA 37
foreach res $resids_CS {
  if {$res <= 37} {
    lappend atomIDs_CS [atomid PROA [expr {$res - 1}] C]
    lappend atomIDs_CS [atomid PROA $res N]
    lappend atomIDs_CS [atomid PROA $res HN]
  } 
  if {$res <= 37} {
    lappend atomIDs_CS [atomid PROB [expr {$res - 1}] C]
    lappend atomIDs_CS [atomid PROB $res N]
    lappend atomIDs_CS [atomid PROB $res HN]
  }
}

foreach atomID [lsort -unique $atomIDs_CS] {
  addatom $atomID
}

proc calcforces {} {
  global initkDC initkCS targetkDC targetkCS targetNumSteps
  global resids_CS atomIDs_CS
  global DCs_exp CSs_exp
  global nu0 sigma11 sigma22 sigma33 sin_delta cos_delta
  global d_xNH_d_xN d_yNH_d_yN d_zNH_d_zN
  global d_xNH_d_xH d_yNH_d_yH d_zNH_d_zH
  global d_xNC_d_xN d_yNC_d_yN d_zNC_d_zN
  global d_xNC_d_xC d_yNC_d_yC d_zNC_d_zC

  set csmax 0
  set csmaxp 0
  set dcmaxp 0
  set dcmax 0

  set kDC 0.
  set kCS 0.
  loadcoords crds
  set timestep [getstep]

  if {$timestep < $targetNumSteps} {
    set scaleFactor [expr {double($timestep) / double($targetNumSteps)}]
    set kDC [expr {$initkDC + $scaleFactor * ($targetkDC - $initkDC)}]
    set kCS [expr {$initkCS + $scaleFactor * ($targetkCS - $initkCS)}]
  } else {
    set kDC $targetkDC
    set kCS $targetkCS
  }

  for {set counter 0} {$counter < [llength $resids_CS]} {incr counter} {
    set atomID_C [lindex $atomIDs_CS [expr {3 * $counter}]]
    set atomID_N [lindex $atomIDs_CS [expr {3 * $counter + 1}]]
    set atomID_H [lindex $atomIDs_CS [expr {3 * $counter + 2}]]
    set DC_exp [lindex $DCs_exp $counter]
    set CS_exp [lindex $CSs_exp $counter]

    # normalized vectors
    set vec_NH [vecnorm [vecsub $crds($atomID_H) $crds($atomID_N)]]
    set vec_NC [vecnorm [vecsub $crds($atomID_C) $crds($atomID_N)]]

    foreach {xNH yNH zNH} $vec_NH {break}
    foreach {xNC yNC zNC} $vec_NC {break}

    # dipolar coupling restraint
    # abs (vec_NH) = 1.
    # eq. s1
    set cos_theta [expr {$zNH}]
    set cos_theta2 [expr {$cos_theta * $cos_theta}]

    # eq. s3
    # alpha = N
    set d_rNH_d_xN [expr {$xNH * $d_xNH_d_xN}]
    set d_rNH_d_yN [expr {$yNH * $d_yNH_d_yN}]
    set d_rNH_d_zN [expr {$zNH * $d_zNH_d_zN}]

    # eq. s2
    # alpha = N
    set d_cosTheta_d_xN [expr {-1. * $zNH * $d_rNH_d_xN}]
    set d_cosTheta_d_yN [expr {-1. * $zNH * $d_rNH_d_yN}]
    set d_cosTheta_d_zN [expr {$d_zNH_d_zN - $zNH * $d_rNH_d_zN}]

    # eq. 4
    # when rigidbond = all, NH distance is kept constant, so
    # the 2nd term of eq. 4 vanishes
    set factor1 [expr {3. * $nu0 * $cos_theta}]

    # alpha = N
    set d_nu_d_xN [expr {$factor1 * $d_cosTheta_d_xN}]
    set d_nu_d_yN [expr {$factor1 * $d_cosTheta_d_yN}]
    set d_nu_d_zN [expr {$factor1 * $d_cosTheta_d_zN}]

    # eq. 1
    set DC_sim [expr {0.5 * $nu0 * (3. * $cos_theta2 - 1.)}]

    set delta_DC 0.
    if {$DC_sim >= 0. && $DC_sim <= $nu0} {
      set delta_DC [expr {$DC_sim - $DC_exp}]
    } elseif {$DC_sim >= [expr {-0.5 * $nu0}] && $DC_sim < 0.} {
      set delta_DC [expr {$DC_sim + $DC_exp}]
    }
    set coeff [expr {-2. * $kDC * $delta_DC}]

    # eq. 4
    set d_nu_d_rN [list $d_nu_d_xN $d_nu_d_yN $d_nu_d_zN]

    # eq. 3
    set d_negativeU_d_rN [vecscale $coeff $d_nu_d_rN]
    set d_negativeU_d_rH [vecscale [expr {-1. * $coeff}] $d_nu_d_rN]
    

	set dcmax [max $dcmax [veclength $d_negativeU_d_rN]]
    if {$dcmax>$dcmaxp} {
		set dcmaxt $atomID_N
		}
	set dcmaxp $dcmax

	set dcmax [max $dcmax [veclength $d_negativeU_d_rH]]
    if {$dcmax>$dcmaxp} {
		set dcmaxt $atomID_H
		}
	set dcmaxp $dcmax
	 

    puts [format "N %4d****%6.2f %6.2f %6.2f, H %4d*****%6.2f %6.2f %6.2f" $atomID_N [lindex $d_negativeU_d_rN 0]   [lindex $d_negativeU_d_rN 1]  [lindex $d_negativeU_d_rN 2] $atomID_H [lindex $d_negativeU_d_rH 0] [lindex $d_negativeU_d_rH 1] [lindex $d_negativeU_d_rH 2]]
    
    addforce $atomID_N $d_negativeU_d_rN
    addforce $atomID_H $d_negativeU_d_rH

    # chemical shift restraint
    # eq. s10
    # e2 = rNC cross_product rNH
    # e2x = yNC * zNH - zNC * yNH
    # e2y = zNC * xNH - xNC * zNH
    # e2z = xNC * yNH - yNC * xNH
    set e2 [vecnorm [veccross $vec_NC $vec_NH]]
    foreach {e2x e2y e2z} $e2 {break}

    # e2 derivatives, i = 2
    # e2x derivatives, e2x = yNC * zNH - zNC * yNH
    # alpha = N
    set d_e2x_d_xN 0.
    set d_e2x_d_yN [expr {$d_yNC_d_yN * $zNH - $zNC * $d_yNH_d_yN}]
    set d_e2x_d_zN [expr {$yNC * $d_zNH_d_zN - $d_zNC_d_zN * $yNH}]

    # alpha = H
    set d_e2x_d_xH 0.
    set d_e2x_d_yH [expr {-1. * $zNC * $d_yNH_d_yH}]
    set d_e2x_d_zH [expr {$yNC * $d_zNH_d_zH}]

    # e2y derivatives, e2y = zNC * xNH - xNC * zNH
    # alpha = N
    set d_e2y_d_xN [expr {$zNC * $d_xNH_d_xN - $d_xNC_d_xN * $zNH}]
    set d_e2y_d_yN 0.
    set d_e2y_d_zN [expr {$d_zNC_d_zN * $xNH - $xNC * $d_zNH_d_zN}]

    # alpha = H
    set d_e2y_d_xH [expr {$zNC * $d_xNH_d_xH}]
    set d_e2y_d_yH 0.
    set d_e2y_d_zH [expr {-1. * $xNC * $d_zNH_d_zH}]

    # e2z derivatives, e2z = xNC * yNH - yNC * xNH
    # alpha = N
    set d_e2z_d_xN [expr {$d_xNC_d_xN * $yNH - $yNC * $d_xNH_d_xN}]
    set d_e2z_d_yN [expr {$xNC * $d_yNH_d_yN - $d_yNC_d_yN * $xNH}]
    set d_e2z_d_zN 0.

    # alpha = H
    set d_e2z_d_xH [expr {-1. * $yNC * $d_xNH_d_xH}]
    set d_e2z_d_yH [expr {$xNC * $d_yNH_d_yH}]
    set d_e2z_d_zH 0.

    # e3
    set oneMinusCosDelta [expr {1. - $cos_delta}]
    set e3 [vecrot $e2 $vec_NH]
    foreach {e3x e3y e3z} $e3 {break}

    # e3 derivatives
    set R11 [expr {$cos_delta + $e2x * $e2x * $oneMinusCosDelta}]
    set R12 [expr {$e2x * $e2y * $oneMinusCosDelta - $e2z * $sin_delta}]
    set R13 [expr {$e2x * $e2z * $oneMinusCosDelta + $e2y * $sin_delta}]
    set R21 [expr {$e2y * $e2x * $oneMinusCosDelta + $e2z * $sin_delta}]
    set R22 [expr {$cos_delta + $e2y * $e2y * $oneMinusCosDelta}]
    set R23 [expr {$e2y * $e2z * $oneMinusCosDelta - $e2x * $sin_delta}]
    set R31 [expr {$e2z * $e2x * $oneMinusCosDelta - $e2y * $sin_delta}]
    set R32 [expr {$e2z * $e2y * $oneMinusCosDelta + $e2x * $sin_delta}]
    set R33 [expr {$cos_delta + $e2z * $e2z * $oneMinusCosDelta}]

    # e3x derivatives
    # e3x = R11 * xNH + R12 * yNH + R13 * zNH
    #     = (cos_delta + e2x ^ 2 * oneMinusCosDelta) * xNH
    #     + (e2x * e2y * oneMinusCosDelta - e2z * sin_delta) * yNH
    #     + (e2x * e2z * oneMinusCosDelta + e2y * sin_delta) * zNH
    # alpha = N
    set d_e3x_d_xN [expr {(2. * $e2x * $d_e2x_d_xN * $oneMinusCosDelta) * $xNH \
                        + $R11 * $d_xNH_d_xN \
                        + ($d_e2x_d_xN * $e2y * $oneMinusCosDelta + $e2x * $d_e2y_d_xN * $oneMinusCosDelta - $d_e2z_d_xN * $sin_delta) * $yNH \
                        + ($d_e2x_d_xN * $e2z * $oneMinusCosDelta + $e2x * $d_e2z_d_xN * $oneMinusCosDelta + $d_e2y_d_xN * $sin_delta) * $zNH \
                        }]
    set d_e3x_d_yN [expr {(2. * $e2x * $d_e2x_d_yN * $oneMinusCosDelta) * $xNH \
                        + ($d_e2x_d_yN * $e2y * $oneMinusCosDelta + $e2x * $d_e2y_d_yN * $oneMinusCosDelta - $d_e2z_d_yN * $sin_delta) * $yNH \
                        + $R12 * $d_yNH_d_yN \
                        + ($d_e2x_d_yN * $e2z * $oneMinusCosDelta + $e2x * $d_e2z_d_yN * $oneMinusCosDelta + $d_e2y_d_yN * $sin_delta) * $zNH \
                        }]
    set d_e3x_d_zN [expr {(2. * $e2x * $d_e2x_d_zN * $oneMinusCosDelta) * $xNH \
                        + ($d_e2x_d_zN * $e2y * $oneMinusCosDelta + $e2x * $d_e2y_d_zN * $oneMinusCosDelta - $d_e2z_d_zN * $sin_delta) * $yNH \
                        + ($d_e2x_d_zN * $e2z * $oneMinusCosDelta + $e2x * $d_e2z_d_zN * $oneMinusCosDelta + $d_e2y_d_zN * $sin_delta) * $zNH \
                        + $R13 * $d_zNH_d_zN}]

    # alpha = H
    set d_e3x_d_xH [expr {(2. * $e2x * $d_e2x_d_xH * $oneMinusCosDelta) * $xNH \
                        + $R11 * $d_xNH_d_xH \
                        + ($d_e2x_d_xH * $e2y * $oneMinusCosDelta + $e2x * $d_e2y_d_xH * $oneMinusCosDelta - $d_e2z_d_xH * $sin_delta) * $yNH \
                        + ($d_e2x_d_xH * $e2z * $oneMinusCosDelta + $e2x * $d_e2z_d_xH * $oneMinusCosDelta + $d_e2y_d_xH * $sin_delta) * $zNH \
                        }]
    set d_e3x_d_yH [expr {(2. * $e2x * $d_e2x_d_yH * $oneMinusCosDelta) * $xNH \
                        + ($d_e2x_d_yH * $e2y * $oneMinusCosDelta + $e2x * $d_e2y_d_yH * $oneMinusCosDelta - $d_e2z_d_yH * $sin_delta) * $yNH \
                        + $R12 * $d_yNH_d_yH \
                        + ($d_e2x_d_yH * $e2z * $oneMinusCosDelta + $e2x * $d_e2z_d_yH * $oneMinusCosDelta + $d_e2y_d_yH * $sin_delta) * $zNH \
                        }]
    set d_e3x_d_zH [expr {(2. * $e2x * $d_e2x_d_zH * $oneMinusCosDelta) * $xNH \
                        + ($d_e2x_d_zH * $e2y * $oneMinusCosDelta + $e2x * $d_e2y_d_zH * $oneMinusCosDelta - $d_e2z_d_zH * $sin_delta) * $yNH \
                        + ($d_e2x_d_zH * $e2z * $oneMinusCosDelta + $e2x * $d_e2z_d_zH * $oneMinusCosDelta + $d_e2y_d_zH * $sin_delta) * $zNH \
                        + $R13 * $d_zNH_d_zH}]
    # e3y derivatives
    # e3y = R21 * xNH + R22 * yNH + R23 * zNH
    #     = (e2y * e2x * oneMinusCosDelta + e2z * sin_delta) * xNH
    #     + (cos_delta + e2y ^ 2 * oneMinusCosDelta) * yNH
    #     + (e2y * e2z * oneMinusCosDelta - e2x * sin_delta) * zNH
    # alpha = N
    set d_e3y_d_xN [expr {($d_e2y_d_xN * $e2x * $oneMinusCosDelta + $e2y * $d_e2x_d_xN * $oneMinusCosDelta + $d_e2z_d_xN * $sin_delta) * $xNH \
                        + $R21 * $d_xNH_d_xN \
                        + (2. * $e2y * $d_e2y_d_xN * $oneMinusCosDelta) * $yNH \
                        + ($d_e2y_d_xN * $e2z * $oneMinusCosDelta + $e2y * $d_e2z_d_xN * $oneMinusCosDelta - $d_e2x_d_xN * $sin_delta) * $zNH \
                        }]
    set d_e3y_d_yN [expr {($d_e2y_d_yN * $e2x * $oneMinusCosDelta + $e2y * $d_e2x_d_yN * $oneMinusCosDelta + $d_e2z_d_yN * $sin_delta) * $xNH \
                        + (2. * $e2y * $d_e2y_d_yN * $oneMinusCosDelta) * $yNH \
                        + $R22 * $d_yNH_d_yN \
                        + ($d_e2y_d_yN * $e2z * $oneMinusCosDelta + $e2y * $d_e2z_d_yN * $oneMinusCosDelta - $d_e2x_d_yN * $sin_delta) * $zNH \
                        }]
    set d_e3y_d_zN [expr {($d_e2y_d_zN * $e2x * $oneMinusCosDelta + $e2y * $d_e2x_d_zN * $oneMinusCosDelta + $d_e2z_d_zN * $sin_delta) * $xNH \
                        + (2. * $e2y * $d_e2y_d_zN * $oneMinusCosDelta) * $yNH \
                        + ($d_e2y_d_zN * $e2z * $oneMinusCosDelta + $e2y * $d_e2z_d_zN * $oneMinusCosDelta - $d_e2x_d_zN * $sin_delta) * $zNH \
                        + $R23 * $d_zNH_d_zN}]

    # alpha = H
    set d_e3y_d_xH [expr {($d_e2y_d_xH * $e2x * $oneMinusCosDelta + $e2y * $d_e2x_d_xH * $oneMinusCosDelta + $d_e2z_d_xH * $sin_delta) * $xNH \
                        + $R21 * $d_xNH_d_xH \
                        + (2. * $e2y * $d_e2y_d_xH * $oneMinusCosDelta) * $yNH \
                        + ($d_e2y_d_xH * $e2z * $oneMinusCosDelta + $e2y * $d_e2z_d_xH * $oneMinusCosDelta - $d_e2x_d_xH * $sin_delta) * $zNH \
                        }]
    set d_e3y_d_yH [expr {($d_e2y_d_yH * $e2x * $oneMinusCosDelta + $e2y * $d_e2x_d_yH * $oneMinusCosDelta + $d_e2z_d_yH * $sin_delta) * $xNH \
                        + (2. * $e2y * $d_e2y_d_yH * $oneMinusCosDelta) * $yNH \
                        + $R22 * $d_yNH_d_yH \
                        + ($d_e2y_d_yH * $e2z * $oneMinusCosDelta + $e2y * $d_e2z_d_yH * $oneMinusCosDelta - $d_e2x_d_yH * $sin_delta) * $zNH \
                        }]
    set d_e3y_d_zH [expr {($d_e2y_d_zH * $e2x * $oneMinusCosDelta + $e2y * $d_e2x_d_zH * $oneMinusCosDelta + $d_e2z_d_zH * $sin_delta) * $xNH \
                        + (2. * $e2y * $d_e2y_d_zH * $oneMinusCosDelta) * $yNH \
                        + ($d_e2y_d_zH * $e2z * $oneMinusCosDelta + $e2y * $d_e2z_d_zH * $oneMinusCosDelta - $d_e2x_d_zH * $sin_delta) * $zNH \
                        + $R23 * $d_zNH_d_zH}]
    # e3z derivatives
    # e3z = R31 * xNH + R32  * yNH + R33 * zNH
    #     = (e2z * e2x * oneMinusCosDelta - e2y * sin_delta) * xNH
    #     + (e2z * e2y * oneMinusCosDelta + e2x * sin_delta) * yNH
    #     + (cos_delta + e2z * e2z * oneMinusCosDelta) * zNH
    # alpha = N
    set d_e3z_d_xN [expr {($d_e2z_d_xN * $e2x * $oneMinusCosDelta + $e2z * $d_e2x_d_xN * $oneMinusCosDelta - $d_e2y_d_xN * $sin_delta) * $xNH \
                        + $R31 * $d_xNH_d_xN \
                        + ($d_e2z_d_xN * $e2y * $oneMinusCosDelta + $e2z * $d_e2y_d_xN * $oneMinusCosDelta + $d_e2x_d_xN * $sin_delta) * $yNH \
                        + (2. * $e2z * $d_e2z_d_xN * $oneMinusCosDelta) * $zNH \
                        }]
    set d_e3z_d_yN [expr {($d_e2z_d_yN * $e2x * $oneMinusCosDelta + $e2z * $d_e2x_d_yN * $oneMinusCosDelta - $d_e2y_d_yN * $sin_delta) * $xNH \
                        + ($d_e2z_d_yN * $e2y * $oneMinusCosDelta + $e2z * $d_e2y_d_yN * $oneMinusCosDelta + $d_e2x_d_yN * $sin_delta) * $yNH \
                        + $R32 * $d_yNH_d_yN \
                        + (2. * $e2z * $d_e2z_d_yN * $oneMinusCosDelta) * $zNH \
                        }]
    set d_e3z_d_zN [expr {($d_e2z_d_zN * $e2x * $oneMinusCosDelta + $e2z * $d_e2x_d_zN * $oneMinusCosDelta - $d_e2y_d_zN * $sin_delta) * $xNH \
                        + ($d_e2z_d_zN * $e2y * $oneMinusCosDelta + $e2z * $d_e2y_d_zN * $oneMinusCosDelta + $d_e2x_d_zN * $sin_delta) * $yNH \
                        + (2. * $e2z * $d_e2z_d_zN * $oneMinusCosDelta) * $zNH \
                        + $R33 * $d_zNH_d_zN}]
    # alpha = H
    set d_e3z_d_xH [expr {($d_e2z_d_xH * $e2x * $oneMinusCosDelta + $e2z * $d_e2x_d_xH * $oneMinusCosDelta - $d_e2y_d_xH * $sin_delta) * $xNH \
                        + $R31 * $d_xNH_d_xH \
                        + ($d_e2z_d_xH * $e2y * $oneMinusCosDelta + $e2z * $d_e2y_d_xH * $oneMinusCosDelta + $d_e2x_d_xH * $sin_delta) * $yNH \
                        + (2. * $e2z * $d_e2z_d_xH * $oneMinusCosDelta) * $zNH \
                        }]
    set d_e3z_d_yH [expr {($d_e2z_d_yH * $e2x * $oneMinusCosDelta + $e2z * $d_e2x_d_yH * $oneMinusCosDelta - $d_e2y_d_yH * $sin_delta) * $xNH \
                        + ($d_e2z_d_yH * $e2y * $oneMinusCosDelta + $e2z * $d_e2y_d_yH * $oneMinusCosDelta + $d_e2x_d_yH * $sin_delta) * $yNH \
                        + $R32 * $d_yNH_d_yH \
                        + (2. * $e2z * $d_e2z_d_yH * $oneMinusCosDelta) * $zNH \
                        }]
    set d_e3z_d_zH [expr {($d_e2z_d_zH * $e2x * $oneMinusCosDelta + $e2z * $d_e2x_d_zH * $oneMinusCosDelta - $d_e2y_d_zH * $sin_delta) * $xNH \
                        + ($d_e2z_d_zH * $e2y * $oneMinusCosDelta + $e2z * $d_e2y_d_zH * $oneMinusCosDelta + $d_e2x_d_zH * $sin_delta) * $yNH \
                        + (2. * $e2z * $d_e2z_d_zH * $oneMinusCosDelta) * $zNH \
                        + $R33 * $d_zNH_d_zH}]

    # e1
    set e1 [vecnorm [veccross $e2 $e3]]
    foreach {e1x e1y e1z} $e1 {break}

    # e1 derivatives, e1 = e2 cross_product e3
    # e1x = e2y * e3z - e2z * e3y
    # alpha = N
    set d_e1x_d_xN [expr {$d_e2y_d_xN * $e3z + $e2y * $d_e3z_d_xN - ($d_e2z_d_xN * $e3y + $e2z * $d_e3y_d_xN)}]
    set d_e1x_d_yN [expr {$d_e2y_d_yN * $e3z + $e2y * $d_e3z_d_yN - ($d_e2z_d_yN * $e3y + $e2z * $d_e3y_d_yN)}]
    set d_e1x_d_zN [expr {$d_e2y_d_zN * $e3z + $e2y * $d_e3z_d_zN - ($d_e2z_d_zN * $e3y + $e2z * $d_e3y_d_zN)}]

    # alpha = H
    set d_e1x_d_xH [expr {$d_e2y_d_xH * $e3z + $e2y * $d_e3z_d_xH - ($d_e2z_d_xH * $e3y + $e2z * $d_e3y_d_xH)}]
    set d_e1x_d_yH [expr {$d_e2y_d_yH * $e3z + $e2y * $d_e3z_d_yH - ($d_e2z_d_yH * $e3y + $e2z * $d_e3y_d_yH)}]
    set d_e1x_d_zH [expr {$d_e2y_d_zH * $e3z + $e2y * $d_e3z_d_zH - ($d_e2z_d_zH * $e3y + $e2z * $d_e3y_d_zH)}]

    # e1y = e2z * e3x - e2x * e3z
    # alpha = N
    set d_e1y_d_xN [expr {$d_e2z_d_xN * $e3x + $e2z * $d_e3x_d_xN - ($d_e2x_d_xN * $e3z + $e2x * $d_e3z_d_xN)}]
    set d_e1y_d_yN [expr {$d_e2z_d_yN * $e3x + $e2z * $d_e3x_d_yN - ($d_e2x_d_yN * $e3z + $e2x * $d_e3z_d_yN)}]
    set d_e1y_d_zN [expr {$d_e2z_d_zN * $e3x + $e2z * $d_e3x_d_zN - ($d_e2x_d_zN * $e3z + $e2x * $d_e3z_d_zN)}]
    
    # alpha = H
    set d_e1y_d_xH [expr {$d_e2z_d_xH * $e3x + $e2z * $d_e3x_d_xH - ($d_e2x_d_xH * $e3z + $e2x * $d_e3z_d_xH)}]
    set d_e1y_d_yH [expr {$d_e2z_d_yH * $e3x + $e2z * $d_e3x_d_yH - ($d_e2x_d_yH * $e3z + $e2x * $d_e3z_d_yH)}]
    set d_e1y_d_zH [expr {$d_e2z_d_zH * $e3x + $e2z * $d_e3x_d_zH - ($d_e2x_d_zH * $e3z + $e2x * $d_e3z_d_zH)}]

    # e1z = e2x * e3y - e2y * e3x
    # alpha = N
    set d_e1z_d_xN [expr {$d_e2x_d_xN * $e3y + $e2x * $d_e3y_d_xN - ($d_e2y_d_xN * $e3x + $e2y * $d_e3x_d_xN)}]
    set d_e1z_d_yN [expr {$d_e2x_d_yN * $e3y + $e2x * $d_e3y_d_yN - ($d_e2y_d_yN * $e3x + $e2y * $d_e3x_d_yN)}]
    set d_e1z_d_zN [expr {$d_e2x_d_zN * $e3y + $e2x * $d_e3y_d_zN - ($d_e2y_d_zN * $e3x + $e2y * $d_e3x_d_zN)}]
    
    # alpha = H
    set d_e1z_d_xH [expr {$d_e2x_d_xH * $e3y + $e2x * $d_e3y_d_xH - ($d_e2y_d_xH * $e3x + $e2y * $d_e3x_d_xH)}]
    set d_e1z_d_yH [expr {$d_e2x_d_yH * $e3y + $e2x * $d_e3y_d_yH - ($d_e2y_d_yH * $e3x + $e2y * $d_e3x_d_yH)}]
    set d_e1z_d_zH [expr {$d_e2x_d_zH * $e3y + $e2x * $d_e3y_d_zH - ($d_e2y_d_zH * $e3x + $e2y * $d_e3x_d_zH)}]

    # eq. s8
    # i = 1, alpha = N
    set d_abs_e1_d_xN [expr {$e1x * $d_e1x_d_xN + $e1y * $d_e1y_d_xN + $e1z * $d_e1z_d_xN}]
    set d_abs_e1_d_yN [expr {$e1x * $d_e1x_d_yN + $e1y * $d_e1y_d_yN + $e1z * $d_e1z_d_yN}]
    set d_abs_e1_d_zN [expr {$e1x * $d_e1x_d_zN + $e1y * $d_e1y_d_zN + $e1z * $d_e1z_d_zN}]

    # i = 1, alpha = H
    set d_abs_e1_d_xH [expr {$e1x * $d_e1x_d_xH + $e1y * $d_e1y_d_xH + $e1z * $d_e1z_d_xH}]
    set d_abs_e1_d_yH [expr {$e1x * $d_e1x_d_yH + $e1y * $d_e1y_d_yH + $e1z * $d_e1z_d_yH}]
    set d_abs_e1_d_zH [expr {$e1x * $d_e1x_d_zH + $e1y * $d_e1y_d_zH + $e1z * $d_e1z_d_zH}]

    # i = 2, alpha = N
    set d_abs_e2_d_xN [expr {$e2x * $d_e2x_d_xN + $e2y * $d_e2y_d_xN + $e2z * $d_e2z_d_xN}]
    set d_abs_e2_d_yN [expr {$e2x * $d_e2x_d_yN + $e2y * $d_e2y_d_yN + $e2z * $d_e2z_d_yN}]
    set d_abs_e2_d_zN [expr {$e2x * $d_e2x_d_zN + $e2y * $d_e2y_d_zN + $e2z * $d_e2z_d_zN}]

    # i = 2, alpha = H
    set d_abs_e2_d_xH [expr {$e2x * $d_e2x_d_xH + $e2y * $d_e2y_d_xH + $e2z * $d_e2z_d_xH}]
    set d_abs_e2_d_yH [expr {$e2x * $d_e2x_d_yH + $e2y * $d_e2y_d_yH + $e2z * $d_e2z_d_yH}]
    set d_abs_e2_d_zH [expr {$e2x * $d_e2x_d_zH + $e2y * $d_e2y_d_zH + $e2z * $d_e2z_d_zH}]

    # i = 3, alpha = N
    set d_abs_e3_d_xN [expr {$e3x * $d_e3x_d_xN + $e3y * $d_e3y_d_xN + $e3z * $d_e3z_d_xN}]
    set d_abs_e3_d_yN [expr {$e3x * $d_e3x_d_yN + $e3y * $d_e3y_d_yN + $e3z * $d_e3z_d_yN}]
    set d_abs_e3_d_zN [expr {$e3x * $d_e3x_d_zN + $e3y * $d_e3y_d_zN + $e3z * $d_e3z_d_zN}]

    # i = 3, alpha = H
    set d_abs_e3_d_xH [expr {$e3x * $d_e3x_d_xH + $e3y * $d_e3y_d_xH + $e3z * $d_e3z_d_xH}]
    set d_abs_e3_d_yH [expr {$e3x * $d_e3x_d_yH + $e3y * $d_e3y_d_yH + $e3z * $d_e3z_d_yH}]
    set d_abs_e3_d_zH [expr {$e3x * $d_e3x_d_zH + $e3y * $d_e3y_d_zH + $e3z * $d_e3z_d_zH}]

    # eq. s7
    # i = 1, alpha = N
    set d_e1z_d_xN [expr {$d_e1z_d_xN - $e1z * $d_abs_e1_d_xN}]
    set d_e1z_d_yN [expr {$d_e1z_d_yN - $e1z * $d_abs_e1_d_yN}]
    set d_e1z_d_zN [expr {$d_e1z_d_zN - $e1z * $d_abs_e1_d_zN}]

    # i = 1, alpha = H
    set d_e1z_d_xH [expr {$d_e1z_d_xH - $e1z * $d_abs_e1_d_xH}]
    set d_e1z_d_yH [expr {$d_e1z_d_yH - $e1z * $d_abs_e1_d_yH}]
    set d_e1z_d_zH [expr {$d_e1z_d_zH - $e1z * $d_abs_e1_d_zH}]

    # i = 2, alpha = N
    set d_e2z_d_xN [expr {$d_e2z_d_xN - $e2z * $d_abs_e2_d_xN}]
    set d_e2z_d_yN [expr {$d_e2z_d_yN - $e2z * $d_abs_e2_d_yN}]
    set d_e2z_d_zN [expr {$d_e2z_d_zN - $e2z * $d_abs_e2_d_zN}]

    # i = 2, alpha = H
    set d_e2z_d_xH [expr {$d_e2z_d_xH - $e2z * $d_abs_e2_d_xH}]
    set d_e2z_d_yH [expr {$d_e2z_d_yH - $e2z * $d_abs_e2_d_yH}]
    set d_e2z_d_zH [expr {$d_e2z_d_zH - $e2z * $d_abs_e2_d_zH}]

    # i = 3, alpha = N
    set d_e3z_d_xN [expr {$d_e3z_d_xN - $e3z * $d_abs_e3_d_xN}]
    set d_e3z_d_yN [expr {$d_e3z_d_yN - $e3z * $d_abs_e3_d_yN}]
    set d_e3z_d_zN [expr {$d_e3z_d_zN - $e3z * $d_abs_e3_d_zN}]

    # i = 3, alpha = H
    set d_e3z_d_xH [expr {$d_e3z_d_xH - $e3z * $d_abs_e3_d_xH}]
    set d_e3z_d_yH [expr {$d_e3z_d_yH - $e3z * $d_abs_e3_d_yH}]
    set d_e3z_d_zH [expr {$d_e3z_d_zH - $e3z * $d_abs_e3_d_zH}]

    # eq. 7
    set CS_sim [expr {$sigma11 * $e1z * $e1z + $sigma22 * $e2z * $e2z + $sigma33 * $e3z * $e3z}]

    # eq. 10
    # alpha = N
    set d_sigma_d_xN [expr {2. * ($sigma11 * $e1z * $d_e1z_d_xN + $sigma22 * $e2z * $d_e2z_d_xN + $sigma33 * $e3z * $d_e3z_d_xN)}]
    set d_sigma_d_yN [expr {2. * ($sigma11 * $e1z * $d_e1z_d_yN + $sigma22 * $e2z * $d_e2z_d_yN + $sigma33 * $e3z * $d_e3z_d_yN)}]
    set d_sigma_d_zN [expr {2. * ($sigma11 * $e1z * $d_e1z_d_zN + $sigma22 * $e2z * $d_e2z_d_zN + $sigma33 * $e3z * $d_e3z_d_zN)}]

    # alpha = H
    set d_sigma_d_xH [expr {2. * ($sigma11 * $e1z * $d_e1z_d_xH + $sigma22 * $e2z * $d_e2z_d_xH + $sigma33 * $e3z * $d_e3z_d_xH)}]
    set d_sigma_d_yH [expr {2. * ($sigma11 * $e1z * $d_e1z_d_yH + $sigma22 * $e2z * $d_e2z_d_yH + $sigma33 * $e3z * $d_e3z_d_yH)}]
    set d_sigma_d_zH [expr {2. * ($sigma11 * $e1z * $d_e1z_d_zH + $sigma22 * $e2z * $d_e2z_d_zH + $sigma33 * $e3z * $d_e3z_d_zH)}]

    set coeff [expr {-2. * $kCS * ($CS_sim - $CS_exp)}]

    set d_sigma_d_rN [list $d_sigma_d_xN $d_sigma_d_yN $d_sigma_d_zN]
    set d_sigma_d_rH [list $d_sigma_d_xH $d_sigma_d_yH $d_sigma_d_zH]
    set d_sigma_d_rC [list [expr {-1. * ($d_sigma_d_xN + $d_sigma_d_xH)}] [expr {-1. * ($d_sigma_d_yN + $d_sigma_d_yH)}] [expr {-1. * ($d_sigma_d_zN + $d_sigma_d_zH)}]]

    set d_negativeU_d_rN [vecscale $coeff $d_sigma_d_rN]
    set d_negativeU_d_rH [vecscale $coeff $d_sigma_d_rH]
    set d_negativeU_d_rC [vecscale $coeff $d_sigma_d_rC]
    

    set csmax [max $csmax [veclength $d_negativeU_d_rC]]
    if {$csmax>$csmaxp} {
		set csmaxt $atomID_C
		}
	set csmaxp $csmax
    set csmax [max $csmax [veclength $d_negativeU_d_rN]]
    if {$csmax>$csmaxp} {
		set csmaxt $atomID_N
		}
	set csmaxp $csmax
    set csmax [max $csmax [veclength $d_negativeU_d_rH]]
    if {$csmax>$csmaxp} {
		set csmaxt $atomID_H
		}
	set csmaxp $csmax
	
	puts [format "C %4d====%6.2f %6.2f %6.2f , N %4d=====%6.2f %6.2f %6.2f, H %4d=====%6.2f %6.2f %6.2f" $atomID_C [lindex $d_negativeU_d_rC 0] [lindex $d_negativeU_d_rC 1] [lindex $d_negativeU_d_rC 2] $atomID_N [lindex $d_negativeU_d_rN 0] [lindex $d_negativeU_d_rN 1] [lindex $d_negativeU_d_rN 2] $atomID_H [lindex $d_negativeU_d_rH 0] [lindex $d_negativeU_d_rH 1] [lindex $d_negativeU_d_rH 2] ]

     addforce $atomID_C $d_negativeU_d_rC
     addforce $atomID_N $d_negativeU_d_rN
     addforce $atomID_H $d_negativeU_d_rH
  }
  puts "now time is $timestep  dcmax is $dcmax at atom $dcmaxt,  csmax is $csmax at atom $csmaxt"
}
