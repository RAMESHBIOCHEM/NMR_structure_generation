
def init():
    """configures C2 symmetry - 
    """

    from atomSel import AtomSel
    #
    # The SymSimulation contains coordinates of copies of the protomer
    # rotated 120 degrees about the z-axis.
    #
    from symSimulation import SymSimulation
    symSim = SymSimulation('symSim',
                           AtomSel("not pseudo"),
                           cloneFirst=False)
    
    from math import pi
    from vec3 import Vec3
    from mat3 import rotVector, Mat3

    N = 2 # homo-oligomer number of units
    
    angle = 2*pi/N  # rotate 360/N degrees about the z-axis

    from utils import char_range
    chars = list( char_range(num=N) )
    for i in range(N):
        symSim.addCopy(rotVector(Vec3(0,0,1), i*angle),
                       Vec3(6,0,0),
                       segidPrefix=chars[i])
        pass

    symSim.numMonomers = N
    symSim.segids = chars
    return symSim
