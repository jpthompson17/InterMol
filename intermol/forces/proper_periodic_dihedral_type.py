import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from abstract_dihedral_type import AbstractDihedralType


class ProperPeriodicDihedralType(AbstractDihedralType):
    """
    Dihedral type with potential given by
    
      V(x) = k * (1 + sign*cos(multiplicity*x - phi)

    where x is the dihedral angle of atoms i,j,k,l in 'protein' convention
    (i.e. x(trans) = 180 degrees)

    NOTE 1: The weight attribute is a weight for 1-4 nonbonded interactions
    between atoms i and l. It doesn't get stored in InterMol's internal
    representation of dihedrals, but it's useful for reading LAMMPS 
    'dihedral_style charmm' parameters.

    """

    __slots__ = ['phi', 'k', 'multiplicity', 'sign', 'weight', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              phi=units.degrees,
                              k=units.kilojoules_per_mole,
                              multiplicity=units.dimensionless,
                              sign=units.dimensionless,
                              weight=units.dimensionless,
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 phi=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 multiplicity=0.0 * units.dimensionless,
                 sign=1 * units.dimensionless,
                 weight=0.0 * units.dimensionless,
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        if (phi != 0 * units.degrees) and (sign != 1 * units.dimensionless):
            # probably didn't mean to set both of these; for now, pass
            pass
        self.phi = phi
        self.k = k
        self.multiplicity = multiplicity
        self.sign = sign
        self.weight = weight


class ProperPeriodicDihedral(ProperPeriodicDihedralType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 phi=0.0 * units.degrees,
                 k=0.0 * units.kilojoules_per_mole,
                 multiplicity=0.0 * units.dimensionless,
                 weight=0.0 * units.dimensionless,
                 improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        ProperPeriodicDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                phi=phi,
                k=k,
                multiplicity=multiplicity,
                weight=weight,
                improper=improper)

    def canonical(self):
        fc = [0.0 * self.k] * 7
        fc[0] = self.k
        fc[self.multiplicity] = self.sign * self.k
    
        return TrigDihedral(self.atom1, self.atom2, self.atom3, self.atom4,
                            self.bondingtype1, self.bondingtype2,
                            self.bondingtype3, self.bondingtype4,
                            self.phi, *fc, self.improper)
