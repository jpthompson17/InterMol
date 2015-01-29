import simtk.unit as units

from intermol.decorators import accepts_compatible_units
from abstract_bond_type import AbstractBondType


class HarmonicBondType(AbstractBondType):
    """
    Bond type with potential given by the function
      
      V(r) = self.k * ( r^order - length^order )^2

    where r is the bond distance.

    NOTE 1: A prefactor may be set at initialization or through a setter
    function. The result of setting the prefactor is that
   
      self.k = prefactor * k

    where k is the force constant specified in __init__. An example of where
    this might be useful is for force fields that specify the harmonic bond
    potential as V = k/2 * (r - length)^2. Note that the prefactor typcically
    should not be queried (hence, "_prefactor").

    """

    __slots__ = ['length', 'k', 'order', 'c', '_prefactor']

    @accepts_compatible_units(None, None, 
                              length=units.nanometers,
                              k=units.kilojoules_per_mole * units.nanometers ** (-2),
                              order=None,
                              c=None,
                              prefactor=None)
    def __init__(self, bondingtype1, bondingtype2, 
                 length=0.0 * units.nanometers,
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 order=1, c=False, prefactor=1):
        AbstractBondType.__init__(self, bondingtype1, bondingtype2, order, c)
        self.length = length
        self.k = prefactor * k
        self._prefactor = prefactor

    @prefactor.setter
    def prefactor(self, prefactor):
        self.k *= prefactor / self._prefactor
        self._prefactor = prefactor


class HarmonicBond(HarmonicBondType):
    """
    stub documentation
    """
    def __init__(self, atom1, atom2, bondingtype1=None, bondingtype2=None, 
                 length=0.0 * units.nanometers,
                 k=0.0 * units.kilojoules_per_mole * units.nanometers ** (-2),
                 order=1, c=False):
        self.atom1 = atom1
        self.atom2 = atom2
        HarmonicBondType.__init__(self, bondingtype1, bondingtype2, 
                length=length,
                k=k,
                order=order, c=c)

    def canonical(self):
        return self
