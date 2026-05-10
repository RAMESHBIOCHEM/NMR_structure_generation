(AB.oligo.refine) 31 % set sel1 [atomselect top "segname 0P1 4P1 CP1 GP1 KP1 OP1 SP1 WP1 and resid 13 to 42"]
atomselect0
>Main< (AB.oligo.refine) 32 % set sel2 [atomselect top "segname 2P1 AP1 EP1 IP1 MP1 QP1 UP1 YP1 and resid 13 to 42"]
atomselect1
>Main< (AB.oligo.refine) 33 % set cent1 [measure center $sel1]
1.7317256927490234 21.128337860107422 14.083964347839355
>Main< (AB.oligo.refine) 34 % set cent2 [measure center $sel2]
-2.13431453704834 -19.99942398071289 -16.28030776977539
>Main< (AB.oligo.refine) 35 % set dist [veclength [vecsub $cent1 $cent2]]
51.26819757243474
>Main< (AB.oligo.refine) 36 % set sel3 [atomselect top "segname 1P1 5P1 DP1 HP1 LP1 PP1 TP1 XP1 and resid 13 to 42"]
atomselect2
>Main< (AB.oligo.refine) 37 % set sel4 [atomselect top "segname 3P1 BP1 FP1 JP1 NP1 RP1 VP1 ZP1 and resid 13 to 42"]
atomselect3
>Main< (AB.oligo.refine) 38 % set cent3 [measure center $cent3]
can't read "cent3": no such variable
>Main< (AB.oligo.refine) 39 % set cent3 [measure center $sel3]
2.34793758392334 17.03179168701172 -21.34060287475586
>Main< (AB.oligo.refine) 40 % set cent4 [measure center $sel4]
1.640423059463501 -16.795202255249023 18.765140533447266
>Main< (AB.oligo.refine) 41 % set dist [veclength [vecsub $cent3 $cent4]]
52.47129453612418
>Main< (AB.oligo.refine) 42 %


set sel1 [atomselect top "protein and (segname 4P1 GP1 OP1 WP1 and resid 10 to 42) or (segname 0P1 CP1 KP1 SP1 and resid 10 to 26) or (segname BP1 JP1 RP1 ZP1 and resid 27 to 42)"]
set sel2 [atomselect top "protein and (segname 2P1 EP1 MP1 UP1 and resid 10 to 42) or (segname AP1 IP1 QP1 YP1 and resid 10 to 26) or (segname 1P1 DP1 LP1 TP1 and resid 27 to 42)"]
set cent1 [measure center $sel1]
set cent2 [measure center $sel2]
set dist [veclength [vecsub $cent1 $cent2]]

step6: 46 distance

###################################################
#
set sel3 [atomselect top "protein and (segname 3P1 FP1 NP1 VP1 and resid 10 to 42) or (segname BP1 JP1 RP1 ZP1 and resid 10 to 26) or (segname AP1 IP1 QP1 YP1 and resid 27 to 42)"]
set sel4 [atomselect top "protein and (segname 5P1 HP1 PP1 XP1 and resid 10 to 42) or (segname 1P1 DP1 LP1 TP1 and resid 10 to 26) or (segname 0P1 CP1 KP1 SP1 and resid 27 to 42)"]
set cent3 [measure center $sel3]
set cent4 [measure center $sel4]
set dist [veclength [vecsub $cent3 $cent4]]

step6: 48.38 distance


##################################################################################################
set sel1 [atomselect top "protein and name CA and (segname 4P1 GP1 OP1 WP1 and resid 10 to 42) or (segname 0P1 CP1 KP1 SP1 and resid 10 to 26) or (segname BP1 JP1 RP1 ZP1 and resid 27 to 42)"]
set sel2 [atomselect top "protein and name CA and (segname 2P1 EP1 MP1 UP1 and resid 10 to 42) or (segname AP1 IP1 QP1 YP1 and resid 10 to 26) or (segname 1P1 DP1 LP1 TP1 and resid 27 to 42)"]
set cent1 [measure center $sel1]
set cent2 [measure center $sel2]
set dist [veclength [vecsub $cent1 $cent2]]

###################################################
#
set sel3 [atomselect top "protein and name CA and (segname 3P1 FP1 NP1 VP1 and resid 10 to 42) or (segname BP1 JP1 RP1 ZP1 and resid 10 to 26) or (segname AP1 IP1 QP1 YP1 and resid 27 to 42)"]
set sel4 [atomselect top "protein and name CA and (segname 5P1 HP1 PP1 XP1 and resid 10 to 42) or (segname 1P1 DP1 LP1 TP1 and resid 10 to 26) or (segname 0P1 CP1 KP1 SP1 and resid 27 to 42)"]
set cent3 [measure center $sel3]
set cent4 [measure center $sel4]
set sel1 [atomselect top "protein and (segname 4P1 GP1 OP1 WP1 and resid 10 to 42) or (segname 0P1 CP1 KP1 SP1 and resid 10 to 26) or (segname BP1 JP1 RP1 ZP1 and resid 27 to 42)"]
set sel2 [atomselect top "protein and (segname 2P1 EP1 MP1 UP1 and resid 10 to 42) or (segname AP1 IP1 QP1 YP1 and resid 10 to 26) or (segname 1P1 DP1 LP1 TP1 and resid 27 to 42)"]
set cent1 [measure center $sel1]
set cent2 [measure center $sel2]
set dist [veclength [vecsub $cent1 $cent2]]set dist [veclength [vecsub $cent3 $cent4]]
