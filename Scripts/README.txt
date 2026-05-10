Generating PSF file and initial coordinates

1) use seq2psf to generate psf files from 3-character sequences:

   generate a protein psf file

        % seq2psf PROT1.seq

          ->  creates PROT1.psf. You can inspect this file with an editor of
              your choice, or

        % idleXplor PROT1.psf

2) use pdb2psf to generate psf file from a pdb file:

     % pdb2psf 1bph.pdb

       -> creates 1bph.psf

    to get all command-line options:
         pdb2psf --help-script

