vmd -dispdev text -e calc_PISEMA.EPRO.final.mod.chainA.tcl

vmd -dispdev text -e calc_PISEMA.EPRO.final.mod.chainB.tcl

paste dc_sim.chainA.xvg dc_sim.chainB.xvg > dc_sim1.dat

paste cs_sim.chainA.xvg cs_sim.chainB.xvg > cs_sim1.dat

awk '{print $1, $2, $4}' dc_sim1.dat > dc_bcsim1.dat

awk '{print $1, $2, $4}' cs_sim1.dat > cs_bcsim1.dat
