>Main< (C2_symmetry) 47 % set com [measure center $all weight mass]
315.4390869140625 313.5528564453125 -1.569042444229126
>Main< (C2_symmetry) 48 % set matrix [transaxis z 180]
{-1.0 -1.2246467991473532e-16 0.0 0.0} {1.2246467991473532e-16 -1.0 0.0 0.0} {0.0 0.0 1.0 0.0} {0.0 0.0 0.0 1.0}
>Main< (C2_symmetry) 49 % $$all moveby [vecscale -1.0 $com]
invalid command name "$atomselect2"
>Main< (C2_symmetry) 50 % $all moveby [vecscale -1.0 $com]
>Main< (C2_symmetry) 51 % $all moveby $com
>Main< (C2_symmetry) 52 % $all writepdb rotate.pdb
Type 'user add key Up {my VMD commands...}' to use this key
>Main< (C2_symmetry) 53 % set all [atomselect top "all"]




>Main< (C2_symmetry) 54 % set com [measure center $all weight mass]
315.4390869140625 313.5528564453125 -1.569042444229126
>Main< (C2_symmetry) 55 % set matrix [transaxis z 90]
{6.123233995736766e-17 -1.0 0.0 0.0} {1.0 6.123233995736766e-17 0.0 0.0} {0.0 0.0 1.0 0.0} {0.0 0.0 0.0 1.0}
>Main< (C2_symmetry) 56 % $all moveby [vecscale -1.0 $com]
>Main< (C2_symmetry) 57 % $all moveby [vecscale -1.0 $com]
>Main< (C2_symmetry) 58 % $all moveby $com
>Main< (C2_symmetry) 59 % $all writepdb rotate2.pdb
>Main< (C2_symmetry) 60 % set all [atomselect top "all"]
atomselect4
>Main< (C2_symmetry) 61 % set com [measure center $all weight mass]
6.107493391027674e-6 2.5656347588665085e-6 -1.6281975945275917e-7
>Main< (C2_symmetry) 62 % set matrix [transaxis z -90]
{6.123233995736766e-17 1.0 0.0 0.0} {-1.0 6.123233995736766e-17 0.0 0.0} {0.0 0.0 1.0 0.0} {0.0 0.0 0.0 1.0}
>Main< (C2_symmetry) 63 % $all move $matrix
>Main< (C2_symmetry) 64 % $all moveby $com
>Main< (C2_symmetry) 65 % $all writepdb rotate.pdb
>Main< (C2_symmetry) 66 % 
