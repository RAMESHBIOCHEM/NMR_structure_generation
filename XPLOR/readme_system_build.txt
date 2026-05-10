1) Calculate structure from random coordinates

  % xplor sdock.py

  (results are in sdock*)

  The protomer structures have suffix .sa, and the full dimer have suffix
  .sa.full.pdb

2) create files with suffix .best for 10 lowest energy structures

  % getBest -symlinks -num 10 sdock_##.sa.stats

  (results are in sdock*.best)

3) The resulting lowest 10 energy structures should have accuracy to the
   reference dimer structure of about 1 angstrom as determined by

targetRMSD -selection 'segid A "" and not name H*' -diffSeq EPRO.pdb `getBest -num 10`

