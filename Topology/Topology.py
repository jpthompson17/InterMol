#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Topology.py

Classes and methods for creating and manipulating molecular system topologies.   

REQUIREMENTS

The SimTK python_units package must be installed See: https://simtk.org/home/python_units




HISTORY

This is based *heavily* on John Chodera's
Pure Python pyopenmm.py code, which is in turn based heavily on Chris Bruns' SimTK PyOpenMM code.

COPYRIGHT

@author Vincent Voelz <vvoelz@gmail.com>
@author John D. Chodera <jchodera@gmail.com>



EXAMPLES


TODO


"""


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import re
import copy
import numpy 

import simtk.unit as units


#=============================================================================================
# EXCEPTIONS
#=============================================================================================

class UnitsException(Exception):
    """
    Exception denoting that an argument has the incorrect units.

    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class ValueException(Exception):
    """
    Exception denoting that an argument has the incorrect value.

    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
    
    

#=============================================================================================
# DECORATORS
#=============================================================================================

#TODO: Do we need to use 'from functools import wraps' to help us here?

def accepts(*types):
    """
    Decorator for class methods that should accept only specified types.

    EXAMPLE

    @accepts(float, int)
    def function(a, b):
        return b*a
    
    """
    def check_accepts(f):
        nargs = (f.func_code.co_argcount - 1) # exclude self
        assert len(types) == nargs, "incorrect number of args supplied in @accepts decorator for class method %s" % (f.func_name)        
        def new_f(*args, **kwds):
            for (a, t) in zip(args[1:], types):
                if a is not None:
                    assert isinstance(a, t), "arg %r does not match %s" % (a,t)
            return f(*args, **kwds)
        new_f.func_name = f.func_name # copy function name
        new_f.func_doc = f.func_doc # copy docstring        
        return new_f
    return check_accepts

def accepts_compatible_units(*units):
    """
    Decorator for class methods that should accept only arguments compatible with specified units.

    Each argument of the function will be matched with an argument of @acceptunits.
    Those arguments of the function that correspond @acceptunits which are not None
    will be checked to ensure they are compatible with the specified units.
    
    EXAMPLE

    @acceptsunits(units.meter, None, units.kilocalories_per_mole)
    def function(a, b, c): pass
    function(1.0 * units.angstrom, 3, 1.0 * units.kilojoules_per_mole)
    
    """
    def check_units(f):
        nargs = (f.func_code.co_argcount - 1) # exclude self
        assert len(units) == nargs, "incorrect number of units supplied in @accepts_compatible_units decorator for class method %s" % (f.func_name)
        def new_f(*args, **kwds):
            for (a, u) in zip(args[1:], units):
                if u is not None:                    
                    assert (a.unit).is_compatible(u), "arg %r does not have units compatible with %s" % (a,u)
            return f(*args, **kwds)
        new_f.func_name = f.func_name # copy function name
        new_f.func_doc = f.func_doc # copy docstring
        return new_f
    return check_units

def returns(rtype):
    """
    Decorator for functions that should only return specific types.
    EXAMPLE

    @returns(int)
    def function(): return 7    
            
    """
    
    def check_returns(f):
        def new_f(*args, **kwds):
            result = f(*args, **kwds)
            assert isinstance(result, rtype), "return value %r does not match %s" % (result,rtype)
            return result
        new_f.func_name = f.func_name # copy function name
        new_f.func_doc = f.func_doc # copy docstring
        return new_f
    return check_returns


    
from Force import *

#=============================================================================================
# Topology base class
#=============================================================================================

class Topology(object):
    """
    This class represents a molecular system topology.  
    
    The Topology object is the outermost container, and stores a list of TopologySystem objects  
    -- The TopologySystem these do most of the heavy lifting.
    """    
    
    def __init__(self, topology=None):
        """
        Create a new Topology object.

        If an Topology object is specified, it will be queried to construct the class.

        """
        self.name = "Untitled"        
        self.molecules      = list()   # molecules[i] is the ith TopologySystem object

    def getName(self):
        """
        Set the name of the Topology.
        
        """
        return self.name

        
    def setName(self, name):
        """
        Set the name of the Topology.
        
        """
        self.name = name
        return 
            

    def getNumMolecules(self):
        """
        Get the number of molecules in the toopology.
        
        """
        return len(self.molecules)
    

    def addMolecule(self, molecule=None):
        """
        Add a molecule (i.e. a TopologySystem object) to the Topology object.
        
        """
        if molecule == None:
            self.molecules.append( TopologySystem() )
        else:
            self.molecules.append( TopologySystem( system=molecule) )
        return


    def delMolecule(self, index):
        """
        Delete a molecule from the Topology object at the specified index.
        
        """
        self.molecules.pop(index)
        return

    def insertMolecule(self, index, molecule=None):
        """
        Insert a molecule into the Topology object at the specified index.
        
        """
        if molecule == None:
            self.molecules.insert( index, TopologySystem() )
        else:
            self.molecules.insert( index, opologySystem( system=molecule) )
        return





class TopologySystem(object):
    """
    A TopologySystem object consists of:
    
        <ol>
        <li>The set of particles in the system</li>
        <li>The forces acting on them</li>
        <li>Pairs of particles whose separation should be constrained to a fixed value</li>
        <li>For periodic systems, the dimensions of the periodic box</li>       
        </ol>
    
        The particles and constraints are defined directly by the System object, while
        forces are defined by objects that extend the Force class.  After creating a
        System, call addParticle() once for each particle, addConstraint() for each constraint,
        and addForce() for each Force.
        
       
    Examples:
    
    >>> system = TopologySystem()

    Add a particle.

    >>> mass = 12.0 * units.amu
    >>> system.addParticle(mass)
    0

    Add a NonbondedForce.

    >>> nonbondedForce = NonbondedForce()
    >>> system.addForce(nonbondedForce)
    0

    Create a deep copy.
    
    >>> import copy
    >>> system_copy = copy.deepcopy(system)


    * A pyopenmm.System-like object stores the atoms in each molecule, their connectivity and forces:

      NOTE that unlike PyOpenMM's System.py, extra metadata information is stored here:
      Atom names, Atom types, along with functions to search through them.


                
    *** USAGE EXAMPLES *** 
    
      
    Topology class should have methods for:
- insert and delete atoms
- insert and delete connections between atoms

- insert and delete atom groups
- merge and split atom groups
 -adding/deleting/joining groups.




- methods for selecting atom groups
   - select by residue
   - select by atomic type
   - select by atom number
   - select by physical location
   
    

    VAV: see  GromacsTopology.py for these:
    - outputting a gromacs top + gro file
    - reading in a topology from a gromacs .itp/.top/.gro
    



- constructing relative free energy topologies given n starting topologies
- building restraints between atoms
- taking a topology and a pdb/gro/etc. and outputting a new topology with the coordinates from the pdb
- creating a 'blank' topology from a pdb or mol2 file (will this be possible?)
- reset arbitrary key pairs in the input information (i.e., we load in the information once from a template file, and then can edit it keyword by keyword as desired).
- setting up A and B states -- that is, a gromacs topology that has a starting state (A) and an ending state (B). (Sort of implied by the relative free energy mention above)
 
    
    
    
    
    
    
    
    *** THINGS THAT WERE ALREADY IMPLEMENTED pyopenmm version*** 
    
    - merging (and splitting?) multiple topologies
    - generating angles from connectivities
    - generating torsions from connectivities
    - generating 1,4's from connectivities
      *  VAV: Force.NonbondedForce.ecreateExceptionsFromBonds()

    VAV:  THIS IS CONTAINED IN THE Force.NonbondedForceParticleInfo class:
    - nonbonded parameters for all atoms  
       - VdW parameters:
             epsilon (float) (kJ/mol)
             sigma (float)
       - partial charge (float)
       - mass (float)
       
    VAV:  THIS IS CONTAINED IN THE MetadataInfo class:
         - involved in free energy calculations (bool)
         - atom name (string)
         - residue name (string)
  
    VAV:  THIS IS CONTAINED IN THE Force.HarmonicBondForceBondInfo class:
     -list of bonds
         - list of parameters: two possibilities -- just a list that can be interpreted by whatever code we want, or exactly what's in the gromacs records
  
    VAV:  THIS IS CONTAINED IN THE Force.HarmonicAngleForceInfo class:
     -list of angles
         - list of parameters: see above
  
    VAV:  THIS IS CONTAINED IN THE Force.PeriodicTorsionForcePeriodicTorsionInfo class:
     -list of torsions
         - list of parameters: see above
   
    VAV:  THIS IS CONTAINED IN THE Force.NonbondedForceExceptionInfo class:
    -list of 1,4's
         - list of parameters: see above
         
    VAV:  THIS IS implemented via the AtomGroup class:
   - list of groups of atoms
      Groups are lists of atom indices, and a name string.
      
    VAV: THIS IS implemented via self.molecules()
      -can represent either single molecules, or systems of molecules

        Q: Can we allow for arbitrary per-atom data to be added later, such as implicit solvent parameters?
           You can either include this information with each atom, or you could think of including all GB data
           as a separate block of information separate from the nonbonded atom terms.I favor the OpenMM strategy
           of having a separate 'Force' object for each force component, so that NonbondedForce terms
           (charges, LJ sigma and epsilon) and metadata on how these nonbonded interactions are computed is separated
           from, say, a generalized Born section which also may have atomic parameters.  The reason for this is that
           the nonbonded parameters will also need associated with them a list of exclusions or exceptions for those
           interactions where the pairwise combining rule needs to be overridden. -John Chodera 5/26/10 7:47 PM 

        VAV sez: For GB parameters, we can use the GB 'Force' objects provided in Force.py.  Otherwise, you can always add new 
        attributes on the fly to the MetadataInfo() objects 


        Lists of "molecules"?
        
        VAV sez:  One requirement was to store a list of "molecules" or objects (since an object may contain multiple molecules
        in the sense of containing multiple disconnected objects); each containing a list of atoms.  In this case, each "molecule"
        should be a separate Topology instance.  The GromacsTopology (and others, potentially) can serve as containers to store
        multiple "molecules"
        
        
        Combination rules to get vdW parameters for unspecified pairs?
        
        MRS sez:  another desired feature is to track nonbonded parameters for specific pairs of atoms to override the above
        JDC sez:  Agreed.  There needs to be the ability to exclude 1,2 and 1,3 pairs and attenuate 1,4 pairs, at the very least.
                  The simplest way is to just include a list of 'exceptions' for nonbonded interactions. -John Chodera 5/26/10 7:49 PM 
        VAV:  This can be implemented via the Force.NonbondedForceExceptionInfo classes and associated methods.


    *** TO DO ***
    
    VAV: A task for MRS's students?
    - select by SMARTS string
     - Hand off the SMARTS string to open eye -- sdf/mol2 file
    

    VAV:  These are NOT implemented yet:  They would require their own derived Topology classes (like GromacsTopology) 
    - reading in a topology from a hetgrpffgen file + pdb (Schrodinger OPLS-AA ?)  VAV: SchrodingerTopology.py ?
    - reading in a topology from Desmond templates(?)  VAV: DesmondTopology.py ?
    - reading in a topology from amber topology files (not necessary yet because of acpypi?) VAV: AmberTopology.py ?

    
    
    *** DISCUSSION - SHOULD THESE TO BE IMPLEMENTED? *** 

    -location: 3d coordinates: triplet of float (nm).
    -velocities?: triplet of float 
    
    - includes the information necessary for the .gro and the .top in gromacs. ???
    
    - list of dummy atoms and parameters for construction of dummy atoms ???



    *** PURPOSELY NOT IMPLEMENTED ***
    
    MRS wrote:  "important note for free energy calculations:Each atom/angle/torsion should be able to have an arbitrary
                number of parameter sets!  For now, it will usually be one or two.  But eventually, we'd like to be able
                to support more.  So it should be a list of lists.
    
    VAV: This can be done for now using multiple objects with different sets of parms.
     
    JDC says:  "A critical question is how alchemical permutations are to be handled by this topology class.
                Should one topology object contain just a single topology, or should there be multiple sets
                of parameters for part or all of the molecule?  I favor just having *one* set of parameters
                defined, and using transformers that take one or two topology objects with identical ordering
                to generate alchemical transformation input files for gromacs, YANK, etc." -John Chodera 5/27/10 6:14 PM 

    Q: Should we potentially include the runfile information, perhaps as a series of string pairs?
       Translators between formats might come later.
    A: I don't think we want runfile info mixed with the topology object. -David Mobley 5/25/10 5:29 PM 

    

    *** TO DO ***
    
    - Very important:  need to lay out what changes need to be made to other files in MM tools
                   to be able to adapt to the new topology definition:
 
        
    """
    
    def __init__(self, system=None):
        """
        Create a new System.

        TODO:
        If a TopologySystem object is specified, its attributes will be queried to construct the class.
        

        """
        
        self.name = ''

        # Set defaults.
        self.masses      = list() # masses[i] is a Quantity denoting the mass of particle i
        self.metadata    = list() # metadata[i] is the MetadataInfo entry for particle i

        self.constraints = list() # constraints[i] is the ith ConstraintInfo entry
        self.forces      = list() # forces[i] is the ith force term
        
        self.atomgroups  = list() # atomgroups[i] is the ith AtomGroup object
 
        
        # VAV: NOT IMPLEMENTED YET
        self.positions   = list() # positions[i] is the ith position 3-tuple (float) 
        self.velocities  = list() # velocities[i] is the ith velocity 3-tuple (float)
        
        
        self.periodicBoxVectors = [ units.Quantity((2.,0.,0.), units.nanometer), units.Quantity((0.,2.,0.), units.nanometer), units.Quantity((0.,0.,2.), units.nanometer) ] # periodic box vectors (only for periodic systems)
        # TODO: Store periodicBoxVectors as units.Quantity(numpy.array([3,3], numpy.float64), units.nanometer)?

        # Populate the system from a provided system, if given.
        if system is not None:
            self._copyDataUsingInterface(self, system)
    
        return

    def _copyDataUsingInterface(self, dest, src):
        """
        Use the public interface to populate 'dest' from 'src'.
        
        """
        #dest.__init__()        
        for index in range(src.getNumParticles()):
            mass = src.getParticleMass(index)            
            dest.addParticle(mass)
        for index in range(src.getNumConstraints()):
            args = src.getConstraintParameters(index)
            dest.addConstraint(*args)
        for index in range(src.getNumForces()):
            force = src.getForce(index)            
            dest.addForce(force)
        box_vectors = src.getPeriodicBoxVectors()
        dest.setPeriodicBoxVectors(*box_vectors)
        
        return
        

    def getName(self):
        """
        Get the name of the TopologySystem.
        """
        return self.name
    
    def setName(self, name):
        """
        Set the name of the TopologySystem.
        """
        self.name = name
        return
    
    
    def getNumParticles(self):
        """
        Get the number of particles in this TopologySystem.
        
        """
        return len(self.masses)

    @accepts_compatible_units(units.amu)
    def addParticle(self, mass):
        """
        Add a particle to the System.

        @param mass   the mass of the particle (in atomic mass units)
        @return the index of the particle that was added

        """

        self.masses.append(mass);
        index = len(self.masses)-1
        self.metadata.append(MetadataInfo(index));  # append a blank metadata descriptor to each particle 

        return index;

    def delParticle(self, index, renumber=True):
        """
        Delete a particle from the System.
        Removes all references (forces, constraints, metadata), bonds) to this particle

        @return the index of the particle that was deleted

        """
        pass

    def insertParticle(self, index, mass, metadata):
        """
        Insert a particle into the System at the specifed insertion index
        Renumbers all references (forces, constraints, metadata), bonds) to this particle accordingly

        @return the index of the particle that was inserted

        """
        pass



    def getParticleMass(self, index):
        """
        Get the mass (in atomic mass units) of a particle.
    
        @param index the index of the particle for which to get the mass

        """        
        return self.masses[index]


    def getParticleMetadata(self, index):
        """
        Get the metadata of a particle.
    
        @param index the index of the particle for which to get the mass

        """        
        return self.metadata[index]



    @accepts_compatible_units(None, units.amu)
    def setParticleMass(self, index, mass):
        """
        Set the mass (in atomic mass units) of a particle.

        @param index  the index of the particle for which to set the mass
        @param mass   the mass of the particle

        """
        masses[index] = mass
        return       
    
    def setParticleMetadata(self, index, atomname=None, atomtype=None, resname=None, resnum=None, \
                                         atomnum=None, atomcharge=None, comment=None, FreeEnergyAtom=False):
        """
        Set the Metadata  of a particle.

        @param index  the index of the particle for which to set the mass

        """
        pass    
    

    def getNumConstraints(self):
        """
        Get the number of distance constraints in this System.

        """
        return len(self.constraints)

    @accepts_compatible_units(None, None, units.nanometer)
    def addConstraint(self, particle1, particle2, distance):
        """
        Add a constraint to the System.
        
        @param particle1 the index of the first particle involved in the constraint
        @param particle2 the index of the second particle involved in the constraint
        @param distance  the required distance between the two particles, measured in nm
        @return the index of the constraint that was added

        """
        if particle1 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")        
        constraint = self.ConstraintInfo(particle1, particle2, distance)
        self.constraints.append(constraint)

        return

    def getConstraintParameters(self, index):
        """
        Get the parameters defining a distance constraint.
        
        @param index     the index of the constraint for which to get parameters
        @return a tuple of (particle1, particle2, distance) for the given constraint index

        """
        constraint = self.constraints[index]
        return (constraint.particle1, constraint.particle2, constraint.distance)

    @accepts_compatible_units(None, None, None, units.nanometer)
    def setConstraintParameters(self, index, particle1, particle2, distance):
        """
        Set the parameters defining a distance constraint.
        
        @param index     the index of the constraint for which to set parameters
        @param particle1 the index of the first particle involved in the constraint
        @param particle2 the index of the second particle involved in the constraint
        @param distance  the required distance between the two particles, measured in nm

        """
        if particle1 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        if particle2 not in range(self.getNumParticles()):
            raise ValueError("particle1 must be in range(0, getNumParticles())")
        constraint = self.ConstraintInfo(particle1,particle2,distance)
        self.constraints[index] = constraint

        return

    def addForce(self, force):
        """
        Add a Force to the System.

        @param force   the Force object to be added
        @return        the index within the System of the Force that was added

        NOTES

        """

        # Append the force.
        self.forces.append(force)
        return len(self.forces)-1

    def getNumForces(self):
        """
        Get the number of Force objects that have been added to the System.

        """
        return len(self.forces)

    def getForce(self, index):
        """
        Get a const reference to one of the Forces in this System.

        @param index  the index of the Force to get

        """
        return self.forces[index]

    def getPeriodicBoxVectors(self):
        """
        Get the vectors which define the axes of the periodic box (measured in nm).  These will affect
        any Force added to the System that uses periodic boundary conditions.

        Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
        x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
        
        @returns a      the vector defining the first edge of the periodic box
        @returns b      the vector defining the second edge of the periodic box
        @returns c      the vector defining the third edge of the periodic box
     
        """        
        return self.periodicBoxVectors

    def setPeriodicBoxVectors(self, a, b, c):
        """
        Set the vectors which define the axes of the periodic box (measured in nm).  These will affect
        any Force added to the System that uses periodic boundary conditions.
        
        Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
        x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
        
        @param a      the vector defining the first edge of the periodic box
        @param b      the vector defining the second edge of the periodic box
        @param c      the vector defining the third edge of the periodic box

        """
        # TODO: Argument checking.
        self.periodicBoxVectors = [a,b,c]
    
    #==========================================================================
    # CONTAINERS
    #==========================================================================

    class ConstraintInfo(object):
        """
        Distance constraint information for particles in System object.

        """

        @accepts_compatible_units(None, None, units.nanometers)
        def __init__(self, particle1, particle2, distance):
            self.particle1 = particle1
            self.particle2 = particle2
            self.distance = distance
            return
        


        
    #==========================================================================
    # PYTHONIC EXTENSIONS
    #==========================================================================    
    
    @property    
    def nparticles(self):
        """
        The number of particles.

        """
        return len(self.masses)

    @property
    def nforces(self):
        """
        The number of force objects in the system.

        """
        return len(self.forces)
        
    @property
    def nconstraints(self):
        """
        The number of interparticle distance constraints defined.

        """
        return len(self.constraints)

    def __str__(self):
        """
        Return an 'informal' human-readable string representation of the System object.

        """
        
        r = ""
        r += "System object\n\n"

        # Show particles.
        r += "Particle masses:\n"
        r += "%8s %24s\n" % ("particle", "mass")
        for index in range(self.getNumParticles()):
            mass = self.getParticleMass(index)
            r += "%8d %24s\n" % (index, str(mass))
        r += "\n"        

        # Show constraints.
        r += "Constraints:\n"
        r += "%8s %8s %16s\n" % ("particle1", "particle2", "distance")
        for index in range(self.getNumConstraints()):
            (particle1, particle2, distance) = self.getConstraintParameters(index)
            r += "%8d %8d %s" % (particle1, particle2, str(distance))
        r += "\n"
        
        # Show forces.
        r += "Forces:\n"
        for force in self.forces:
            r += str(force)
            
        return r

    def __add__(self, other):
        """
        Binary concatenation of two systems.

        The atoms of the second system appear, in order, after the atoms of the first system
        in the new combined system.

        USAGE

        combined_system = system1 + system2

        NOTES

        Both systems must have identical ordering of Force terms.
        Any non-particle settings from the first System override those of the second, if they differ.

        EXAMPLES

        Concatenate two systems in vacuum.

        >>> import testsystems
        >>> [system1, coordinates1] = testsystems.LennardJonesFluid()
        >>> system1 = System(system1) # convert to pure Python system
        >>> [system2, coordinates2] = testsystems.LennardJonesFluid()
        >>> system2 = System(system2) # convert to pure Python system        
        >>> combined_system = system1 + system2

        """
        system = copy.deepcopy(self)
        system += other
        return system
        
    def __iadd__(self, other):
        """
        Append specified system.

        USAGE

        system += additional_system

        NOTES

        Both systems must have identical ordering of Force terms.
        Any non-particle settings from the first System override those of the second, if they differ.
        
        EXAMPLES

        Append atoms from a second system to the first.
    
        >>> import testsystems
        >>> [system1, coordinates1] = testsystems.LennardJonesFluid()
        >>> system1 = System(system1) # convert to pure Python system        
        >>> [system2, coordinates2] = testsystems.LennardJonesFluid()
        >>> system2 = System(system2) # convert to pure Python system        
        >>> system1 += system2

        """
        # Check to make sure both systems are compatible.
        if not isinstance(other, System):
            raise ValueError("both arguments must be System objects")
        if (self.nforces != other.nforces):
            raise ValueError("both System objects must have identical number of Force classes")
        for (force1,force2) in zip(self.forces, other.forces):
            if type(force1) != type(force2):
                raise ValueError("both System objects must have identical ordering of Force classes")

        # TODO: Make sure other system is Pythonic instance.

        # Combine systems.
        for mass in other.masses:
            mass_copy = copy.deepcopy(mass)
            self.masses.append(mass_copy)
        for constraint in other.constraints:
            constraint_copy = copy.deepcopy(constraint)
            self.constraints.append(constraint_copy)
        for (force1, force2) in zip(self.forces, other.forces):
            force2_copy = copy.deepcopy(force2)
            offset = self.nparticles
            force1._appendForce(force2_copy, offset)
                
        return

    #==========================================================================
    # Methods for selecting atom groups based on Metadata
    #==========================================================================
   
    def selectAtomGroup(self, atomname=None, atomtype=None, resname=None, resnum=None, \
                              atomnum=None, atomcharge=None, freeEnergyAtom=False):
        
        """Select particles according to specified metadata attributes.  Selections can be lists,
        or individual entries.
        
        Examples:
        
        Select all atomnames that are either N, CA or C, and in residue 4

        >>> selection = selectAtomGroup(atomname=['N','CA','C'], resnum=4)
        >>> selection = selectAtomGroup(atomname=['N','CA','C'], resnum=[4])
        
        RETURNS an AtomGroup object
        """
       
        # Convert the specification to lists (or None)
        if not isinstance(atomname,list) and (atomname!=None):
            AtomNameList = [ atomname ]
        if not isinstance(atomtype,list) and (atomtype!=None):
            AtomTypeList = [ atomtype ]
        if not isinstance(resname,list) and (resname!=None):
            ResidueNameList = [ resname ]
        if not isinstance(resnum,list) and (resnum!=None):
            ResidueNumList = [ resnum ]
        if not isinstance(atomcharge,list) and (atomcharge!=None):
            AtomChargeList = [ atomcharge ]
        if not isinstance(freeEnergyAtom,list) and (freeEnergyAtom!=None):
            freeEnergyList = [ freeEnergyAtom ]
        
        group = AtomGroup()
        for metadatum in self.metadata:
            Keep = True
            Keep *= ((metadatum.atomname in AtomNameList) or AtomNameList==None)
            Keep *= ((metadatum.atomtype in AtomTypeList) or AtomTypeList==None)
            Keep *= ((metadatum.resname in ResidueNameList) or ResidueNameList==None)
            Keep *= ((metadatum.resnum in ResidueNumList) or ResidueNumList==None)
            Keep *= ((metadatum.atomcharge in AtomChargeList) or AtomChargeList==None)
            Keep *= ((metadatum.freeEnergyAtom in freeEnergyList) or freeEnergyList==None)
            if Keep:
                group.addIndex(metadatum.particle)
            
        
    #==========================================================================
    # Methods for get and set of MetaData quantities
    #==========================================================================
    
    # by particle (index)
    
    def getAtomNameByParticle(self):
        pass

    def setAtomNameByParticle(self):
        pass

    def getAtomTypeByParticle(self):
        pass

    def setAtomTypeByParticle(self):
        pass
    
    def getResidueNameByParticle(self):
        pass

    def setResidueNameByParticle(self):
        pass

    def getResidueNumByParticle(self):
        pass

    def setResidueNumByParticle(self):
        pass
    
    def getAtomChargeByParticle(self):
        pass

    def setAtomChargeByParticle(self):
        pass
    
    def getCommentByParticle(self):
        pass

    def setCommentByParticle(self):
        pass

    def getFreeEnergyAtomByParticle(self):
        pass

    def setFreeEnergyAtomByParticle(self):
        pass
    
    def getAtomNameByParticle(self):
        pass
    
    
    # by atomgroup (index)

    def setAtomNameByAtomGroup(self):
        pass

    def getAtomTypeByAtomGroup(self):
        pass

    def setAtomTypeByAtomGroup(self):
        pass
    
    def getResidueNameByAtomGroup(self):
        pass

    def setResidueNameByAtomGroup(self):
        pass

    def getResidueNumByAtomGroup(self):
        pass

    def setResidueNumByAtomGroup(self):
        pass
    
    def getAtomChargeByAtomGroup(self):
        pass

    def setAtomChargeByAtomGroup(self):
        pass
    
    def getCommentByAtomGroup(self):
        pass

    def setCommentByAtomGroup(self):
        pass

    def getFreeEnergyAtomByAtomGroup(self):
        pass

    def setFreeEnergyAtomByAtomGroup(self):
        pass

    
    
class MetadataInfo(object):
    """
    A container class for storing atom names, types, residues
    """

    def __init__(self, particle, atomname=None, atomtype=None, resname=None, resnum=None, \
                                  atomnum=None, atomcharge=None, comment=None, freeEnergyAtom=False):
        """
        ;     nr type        resnr residue   atom   cgnr charge     mass       typeB chargeB massB
               1 amber99_39      1   NPRO      N      1 -0.202     14.01      ; qtot -0.202
        """
        self.particle = particle
        self.atomname = atomname
        self.atomtype = atomtype
        self.resname = resname
        self.resnum = resnum
        self.atomnum = atomnum
        self.atomcharge = atomcharge
        self.comment = comment
        self.freeEnergyAtom = freeEnergyAtom  # boolean
        return

    def getAtomName(self):
        pass

    def setAtomName(self):
        pass
    
    def getAtomType(self):
        pass

    def setAtomType(self):
        pass
    
    def getResidueName(self):
        pass

    def setResidueName(self):
        pass

    def getResidueNum(self):
        pass

    def setResidueNum(self):
        pass
    
    def getAtomCharge(self):
        pass

    def setAtomCharge(self):
        pass
    
    def getComment(self):
        pass

    def setComment(self):
        pass

    def getFreeEnergyAtom(self):
        pass

    def setFreeEnergyAtom(self):
        pass

    
    
    
class AtomGroup(object):
    """
    A class that stores a list of atom indices and an atomgroup name.
    Selection methods will return AtomGroup objects as the result of a selection query.
    
    Metadata properties can also be set using AtomGroup objects 
    
    """
    
    def __init__(self, atomgroup=None):
        """Initialize an AtomGroup class.  If an input class is provided, instantiate using the information from it."""
        
        self.indices = list()
        self.name = None
        
        if atomgroup:
            self.indices = copy.deepcopy(atomgroup.indices)
            self.name = copy.deepcopy(atomgroup.name)
          
        return

    def addIndex(self, index):
        self.indices.append(index)



    