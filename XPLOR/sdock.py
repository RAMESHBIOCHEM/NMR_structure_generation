
xplor.requireVersion("3.4")

#
# strict symmetry docking/structure calculation of homotrimer
#

# this checks for typos on the command-line. User-customized arguments can
# also be specified.
#
(opts,args) = xplor.parseArguments(["quick"])

quick=False
for opt in opts:
    if opt[0]=="quick":  #specify -quick to just test that the script runs
        quick=True
        pass
    pass


# filename for output structures. This string must contain the STRUCTURE
# literal so that each calculated structure has a unique name. The SCRIPT
# literal is replaced by this filename (or stdin if redirected using <),
# but it is optional.
#
outFilename = "SCRIPT_STRUCTURE.sa"
numberOfStructures=500  #usually you want to create at least 20 

if quick:
    numberOfStructures=2
    pass

# protocol module has many high-level helper functions.
#
import protocol

protocol.initRandomSeed()   #set random seed - by time

command = xplor.command

# generate PSF data from sequence and initialize the correct parameters.
#
from psfGen import seqToPSF
seqToPSF('EPRO.seq',startResid=12)

# generate random extended initial structure with correct covalent geometry
#
protocol.genExtendedStructure()

#all symmetry-related subunits contained in symSim
# symSimSetup module contains the symmetry operations for the current
# calculation.
import symSimSetup
symSim = symSimSetup.init()  

numMonomers = symSim.numMonomers 
segids = symSim.segids # segids of the subunits defined in the SymSimulation

#
# a PotList contains a list of potential terms. This is used to specify which
# terms are active during refinement.
#
from potList import PotList
potList = PotList()

# parameters to ramp up during the simulated annealing protocol
#
from simulationTools import MultRamp, StaticRamp, InitialParams, IVMAction

rampedParams=[]
highTempParams=[]

# IVM setup
#   the IVM is used for performing dynamics and minimization in torsion-angle
#   space, and in Cartesian space.
#   In strict symmetry calculations, only the protomer coordinates ( in
#   xplor.simulation) are manipulated by the IVM objects
#
from ivm import IVM
dyn  = IVM()
minc = IVM() # minc used for final cartesian minimization

# initialize ivm topology for torsion-angle dynamics


#
# 
#


# here, we instantiate energy terms involving the full construct
# the default set of coordinates is now those of the mull multimer.
#
import simulation
simulation.makeCurrent( symSim )

# compare atomic Cartesian rmsd with a reference structure
#  backbone and heavy atom RMSDs will be printed in the output
#  structure files
#
#from posDiffPotTools import create_PosDiffPot
#refRMSD = create_PosDiffPot("refRMSD",
#                            AtomSel("name CA",symSim),
#                            pdbFile='EPRO.pdb',
#                            cmpSel="not name H*")
#

# set up NOE potential for distance restraints
# restaints here were identified as intermolecular, and thus
# have segid specified, For restraints which may be inter- or intramolecular
# the segid should be omitted, and numMono should be set to numMonomers.
#
noeInter=PotList('noeInter')
potList.append(noeInter)
from noePotTools import create_NOEPot
for (name,scale,file) in [('all',1,"inter.tbl"),
                          #add entries for additional tables
                          ]:
    pot = create_NOEPot(name,file)
    # pot.setPotType("soft") # if you think there may be bad NOEs
    pot.setScale(scale)
    noeInter.append(pot)
rampedParams.append( MultRamp(2,30, "noeInter.setScale( VALUE )") )




# gyration volume term - pack the subunits together - this assumes
# that the overall structure is fairly compact. If not, a subset of atoms
# might be specified which are thought to have this property.
#
from gyrPotTools import create_GyrPot
gyr = create_GyrPot("Vgyr",
#                    "resid A:B" # selection should exclude disordered tails
                    ) 
potList.append(gyr)
rampedParams.append( MultRamp(.002,1,"gyr.setScale(VALUE)") )


# HBPot - knowledge-based hydrogen bond term
# this acts on both inter- and intramolecular hydrogen bonds
#
# maybe best to use during refinement only
from hbPotTools import create_HBPot
hb = create_HBPot('hb')
hb.setScale(2.5)
potList.append( hb )

#
# setup parameters for atom-atom repulsive term. (van der Waals-like term)
# Due to symmetry, only a subset of subunits need be specified in the selPairs
# argument.
#
from repelPotTools import create_RepelPot,initRepel
#sufficient to compute interactions between one subunit and each of the others
repel = create_RepelPot('repel',
                        selPairs=[('segid A',f'segid {segid}') for segid
                                  in segids])
potList.append(repel)
rampedParams.append( StaticRamp("initRepel(repel,use14=False)") )
rampedParams.append( MultRamp(.004,4,  "repel.setScale( VALUE*numMonomers)") )
# nonbonded interaction only between CA atoms
highTempParams.append( StaticRamp("""initRepel(repel,
                                               use14=True,
                                               scale=0.004,
                                               repel=1.2,
                                               moveTol=45,
                                               interactingAtoms='name CA'
                                               )""") )



# the following terms act only on the protomer coordinates
simulation.makeCurrent( xplor.simulation )

# set up intra-subunit distance restraints
# segids should not be specified in the restraints, as the protomer atoms
# have none - none were specified in seqToPSF
noeIntra=PotList('noeIntra')
potList.append(noeIntra)
from noePotTools import create_NOEPot
for (name,scale,file) in [('hb',1,"hbond.tbl"),
                          #add entries for additional tables
                          ]:
    pot = create_NOEPot(name,file)
    # pot.setPotType("soft") # if you think there may be bad NOEs
    pot.setScale(scale)
    noeIntra.append(pot)
rampedParams.append( MultRamp(2,30, "noeIntra.setScale( VALUE*numMonomers )") )


# Set up dihedral angle restraints 
dihedralRestraintFilename="dihed.tbl"
from dihedralPotTools import create_DihedralPot
dihePot = create_DihedralPot('dihePot',dihedralRestraintFilename)
potList.append( dihePot )
highTempParams.append( StaticRamp("dihePot.setScale(10*numMonomers)") )
rampedParams.append( StaticRamp("dihePot.setScale(200*numMonomers)") )

#Torsion angle knowledge-based database potential
#
from torsionDBPotTools import create_TorsionDBPot
torsionDB = create_TorsionDBPot('torsionDB', system='protein')
potList.append( torsionDB )
rampedParams.append( MultRamp(.002,2,"torsionDB.setScale(VALUE*numMonomers)") )

# Selected 1-4 nonbonded interactions to avoid eclipsed conformations of
# protons on sidechain termini
import torsionDBPotTools
repel14 = torsionDBPotTools.create_Terminal14Pot('repel14')
potList.append(repel14)
highTempParams.append(StaticRamp("repel14.setScale(0)"))
rampedParams.append(MultRamp(0.004, 4, "repel14.setScale(VALUE*numMonomers)"))

# convalent energy terms acting on the protomer
from xplorPot import XplorPot
potList.append( XplorPot("BOND") )
potList['BOND'].setScale( numMonomers )
potList.append( XplorPot("ANGL") )
potList['ANGL'].setThreshold( 5 )
rampedParams.append( MultRamp(0.4,1,
                              "potList['ANGL'].setScale(VALUE*numMonomers)") )
potList.append( XplorPot("IMPR") )
potList['IMPR'].setThreshold( 5 )
rampedParams.append( MultRamp(0.1,1,
                              "potList['IMPR'].setScale(VALUE*numMonomers)") )
      


# Give atoms uniform weights, except for the anisotropy axis
#
protocol.massSetup()


#configure dyn IVM object to act in torsion angle space
protocol.torsionTopology(dyn)

#configure minc IVM object to act in Cartesian space
protocol.cartesianTopology(minc)



# object which performs simulated annealing
#
from simulationTools import AnnealIVM
init_t  = 3500.     # Need high temp and slow annealing to converge
cool = AnnealIVM(initTemp =init_t,
                 finalTemp=25,
                 tempStep =12.5,
                 ivm=dyn,
                 rampedParams = rampedParams)

#cart_cool is for optional cartesian-space cooling
cart_cool = AnnealIVM(initTemp =init_t,
	              finalTemp=25,
		      tempStep =12.5,
                      ivm=minc,
                      rampedParams = rampedParams)

def calcOneStructure(loopInfo):
    """ this function calculates a single structure, performs analysis on the
    structure, and then writes out a pdb file, with remarks.
    """

    atomPos0 = xplor.simulation.atomPosArr()
    atomPosMin = xplor.simulation.atomPosArr()
    energyMin =  1e30	# big number
    k=0		# set up initial minimaizing loop

    # this loop attempts 5 docking minimnizations - from different, randomized
    # subunit positions/orientations, and with different randomized internal
    # torsion angles.
    while k < 5:
        # reset coordinates to original values.
        xplor.simulation.setAtomPosArr(atomPos0)
        
        #randomize internal torsion angles
        from monteCarlo import randomizeTorsions
        randomizeTorsions(dyn)

        #randomize the protomer orientation and position
        # note that symmetry of the full set of coordinates in symSim is
        # always maintained.        
        from atomAction import randomizeDomainPos
        randomizeDomainPos("not pseudo")
        
        
        # set torsion angles from restraints
        #
        from torsionTools import setTorsionsFromTable
        setTorsionsFromTable(dihedralRestraintFilename)

        #initially use only distance restraints and nonbonded repulsion
        #in an energy minimization.
        protocol.initMinimize(dyn,potList=[noeInter,repel],
                              printInterval=50)
        dyn.run() ; dyn.run()
        
        # next, do minimization using all energy terms.
        #
        protocol.initMinimize(dyn,potList=potList,
                              printInterval=50)
        dyn.run() ; dyn.run()
        
        dyn.run() ; dyn.run()

        print( 'iter: %d    energy: %f' % (k,potList.calcEnergy()))

        #save this set of protomer coords if it produces the lowest energy.
        if potList.calcEnergy() < energyMin :
            atomPosMin = xplor.simulation.atomPosArr()
            energyMin = potList.calcEnergy()
            pass
        
        k += 1
        pass

    #use the lowest energy pose
    xplor.simulation.setAtomPosArr(atomPosMin)
    protocol.fixupCovalentGeom(maxIters=100,useVDW=1)

    #the .init file contains coordinates from the initial docking calculation
    #with fixed covalent geometry.
    protocol.writePDB(loopInfo.filename()+".init")

    # initialize parameters for high temp dynamics.
    InitialParams( rampedParams )
    # high-temp dynamics setup - only need to specify parameters which
    #   differfrom initial values in rampedParams
    InitialParams( highTempParams )

    # high temp dynamics
    #
    protocol.initDynamics(dyn,
                          potList=potList, # potential terms to use
                          bathTemp=init_t,
                          initVelocities=1,
                          finalTime=100,   # stops at 100ps or 1000 steps
                          numSteps=1000,   # whichever comes first
                          printInterval=100)

    dyn.setETolerance( init_t/100 )  #used to det. stepsize. default: t/1000 
    dyn.run()

    # initialize integrator for simulated annealing
    #
    protocol.initDynamics(dyn,
                          numSteps=100,       #at each temp: 100 steps or
                          finalTime=.2 ,       # .2ps, whichever is less
                          printInterval=100)

    # perform simulated annealing
    #
    cool.run()
              
              
    # final torsion angle minimization
    #
    protocol.initMinimize(dyn,
                          printInterval=50)
    dyn.run()

    # optional cooling in Cartesian coordinates
    #
    protocol.initDynamics(minc,
                          potList=potList,
                          numSteps=100,       #at each temp: 100 steps or
                          finalTime=.4 ,       # .2ps, whichever is less
                          printInterval=100)
    #cart_cool.run()
    # final all- atom minimization
    #
    protocol.initMinimize(minc,
                          potList=potList,
                          dEPred=10)
    minc.run()

    #do analysis and write structure when function returns
    # the coordinates of the full multimer are written to a file with suffix
    # .full
    protocol.writePDB(loopInfo.filename()+".full.pdb",
                      selection=AtomSel("all",symSim))
    pass


from simulationTools import StructureLoop, FinalParams
StructureLoop(numStructures=numberOfStructures,
              pdbTemplate=outFilename,
              structLoopAction=calcOneStructure,
              doWriteStructures=True,
              genViolationStats=True,
              averageTopFraction=0.5, #report stats on best 50% of structs
              averageContext=FinalParams(rampedParams),
              #averageCrossTerms=refRMSD,
              averageFitSel=None, #do not fit subunit coords (this would fit
                                  #protomer coorda, and thus change packing).
              averagePotList=potList).run()
