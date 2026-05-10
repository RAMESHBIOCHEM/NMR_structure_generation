vmd -dispdev text -e dihedral_amber_file.tcl

#grep -v "DIHEDRAL  ZFC 0.0" dihe.txt > temp.txt

#grep -v "IMPROPER   ZFC 120.0" temp.txt > dihe.txt

#rm temp.txt

sed -i -e 's/ZFC/$FC/g' Dihedral.ref


#vmd -dispdev text -e rest_pro_lipid_amber.tcl

#vmd -dispdev text -e rest_pro_sc_bb.tcl
