import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from abstract_dihedral_type import AbstractDihedralType


class RbDihedralType(AbstractDihedralType):
    """
    Dihedral type with potential given by the Ryckaert-Bellemans (RB) function:
    
      V_rb(x) = sum( C_n * (cos(x))^n, n = 0,...,6 )                    [1]
    
    For the 'polymer' convention (default),
    
      x = psi = phi - 180, with psi(trans) = 0
    
    For the 'protein' convention,
    
      x = phi, with phi(trans) = 180
    
    Conversion between conventions can be accomplished by multiplying each
    coefficient C_n by (-1)^n and replacing psi with phi (or vice versa)
    in Eq. 1 above. (See the GROMACS 5.0.2 Manual, Section 4.2.13.)
    
    Examples: The polymer convention is used for RB (func. type 3) dihedrals
    in GROMACS. The protein convention is used for the 'multi/harmonic' 
    dihedral style in LAMMPS.

    NOTE 1: The convention set on construction of a RbDihedralType (again,
    'polymer' by default) should be consistent with its initial set of 
    C* coefficients. Subsequent changes made to the convention will 
    automatically flip the sign of C1, C3, and C5.

    """

    __slots__ = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', '_convention', 'improper']

    @accepts_compatible_units(None, None, None, None, 
                              C0=units.kilojoules_per_mole,
                              C1=units.kilojoules_per_mole,
                              C2=units.kilojoules_per_mole,
                              C3=units.kilojoules_per_mole,
                              C4=units.kilojoules_per_mole,
                              C5=units.kilojoules_per_mole,
                              C6=units.kilojoules_per_mole,
                              convention=None, 
                              improper=None)
    def __init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                 C0=0.0 * units.kilojoules_per_mole,
                 C1=0.0 * units.kilojoules_per_mole,
                 C2=0.0 * units.kilojoules_per_mole,
                 C3=0.0 * units.kilojoules_per_mole,
                 C4=0.0 * units.kilojoules_per_mole,
                 C5=0.0 * units.kilojoules_per_mole,
                 C6=0.0 * units.kilojoules_per_mole,
                 convention='polymer',
                 improper=False):
        AbstractDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, improper)
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        self.C5 = C5
        self.C6 = C6
        self._convention = convention

    @property
    def convention(self):
        return self._convention

    @convention.setter
    def convention(self, convention):
        if not (convention == 'polymer' or convention == 'protein'):
            e = ValueError("Invalid convention for RbDihedralType: {0}".format(convention))
            logger.exception(e)

        # No action if convention is the same
        if (convention == self._convention):
            return

        # Otherwise, switch conventions
        self._convention = convention
        self.C1 *= -1
        self.C3 *= -1
        self.C5 *= -1

class RbDihedral(RbDihedralType):
    """
    Dihedral interaction of Ryckaert-Bellemans type (see comments in 
    RbDihedralType for more information).
    """
    def __init__(self, atom1, atom2, atom3, atom4, bondingtype1=None, bondingtype2=None, bondingtype3=None, bondingtype4=None, 
                 C0=0.0 * units.kilojoules_per_mole,
                 C1=0.0 * units.kilojoules_per_mole,
                 C2=0.0 * units.kilojoules_per_mole,
                 C3=0.0 * units.kilojoules_per_mole,
                 C4=0.0 * units.kilojoules_per_mole,
                 C5=0.0 * units.kilojoules_per_mole,
                 C6=0.0 * units.kilojoules_per_mole,
                 convention='polymer',
                 improper=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        RbDihedralType.__init__(self, bondingtype1, bondingtype2, bondingtype3, bondingtype4, 
                C0=C0,
                C1=C1,
                C2=C2,
                C3=C3,
                C4=C4,
                C5=C5,
                C6=C6,
                convention=convention,
                improper=improper)

    def canonical(self):
        """
        Return a TrigDihedral equivalent to this RbDihedral.
        """

        # Convert to same convention as TrigDihedral (protein/phi convention)
        init_convention = self.convention
        self.convention = 'protein'

        # Calculate TrigDihedral parameters
        phi_phase = 0 * units.degrees

        fc0 = self.C0 + 0.5*self.C2 + 0.375*self.C4 + 0.3125 *self.C6
        fc2 =           0.5*self.C2 + 0.5  *self.C4 + 0.46875*self.C6
        fc4 =                         0.125*self.C4 + 0.1875 *self.C6
        fc6 =                                         0.03125*self.C6
        
        fc1 = self.C1 + 0.75*self.C3 + 0.625 *self.C5
        fc3 =           0.25*self.C3 + 0.3125*self.C5
        fc5 =                          0.0625*self.C5
        
        # Restore to initial convention
        self.convention = init_convention

        return TrigDihedral(self.atom1, self.atom2, self.atom3, self.atom4,
                            self.bondingtype1, self.bondingtype2,
                            self.bondingtype3, self.bondingtype4,
                            phi_phase, fc0, fc1, fc2, fc3, fc4, fc5, fc6,
                            self.improper)
