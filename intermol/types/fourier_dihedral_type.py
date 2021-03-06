import sys
sys.path.append('..')

from intermol.decorators import *
from abstract_dihedral_type import *

class FourierDihedralType(AbstractDihedralType):

    @accepts_compatible_units(None,
            None,
            None,
            None,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole,
            units.kilojoules_per_mole)

    def __init__(self, atom1, atom2, atom3, atom4, c1, c2, c3, c4):
        """
        """
        AbstractDihedralType.__init__(self, atom1, atom2, atom3, atom4)
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4

