proc vec_Epro {} {

## tilt angle check
	set out [open "vec_Edimer.step8.1.dat" w]
	
	for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
	
	
                set a1 [atomselect top "protein and resid 25 to 37 and segname PROA and name CA" frame $i]
                set a2 [atomselect top "protein and resid 12 to 25 and segname PROA and name CA" frame $i]

                set cm1 [measure center $a1 weight mass]
                set cm2 [measure center $a2 weight mass]
                set vec [vecsub $cm1 $cm2]

                set vz [list 0 0 1]

                set dp [vecdot $vz $vec]
                set cos [expr $dp/[veclength $vec] ]
                set proa [expr acos($cos)*(180/3.14)]
                  
                ##################################################################3
		set b1 [atomselect top "protein and resid 25 to 37 and segname PROB and name CA" frame $i]
                set b2 [atomselect top "protein and resid 12 to 25 and segname PROB and name CA" frame $i]

                set cmb1 [measure center $b1 weight mass]
                set cmb2 [measure center $b2 weight mass]
                set vecb [vecsub $cmb1 $cmb2]

                set vz [list 0 0 1]

                set dpb [vecdot $vz $vecb]
                set cos [expr $dpb/[veclength $vecb] ]

                set prob [expr acos($cos)*(180/3.14)]

		
		puts $out "$i $proa $prob"
		
		$a1 delete
		$a2 delete
		$b1 delete
		$b2 delete
	}
	
	close $out
}  
