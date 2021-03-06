from intermol.decorators import *
from abstract_pair import *

class LJ2PairCR1(AbstractPair):
    __slots__ = ['fudgeQQ', 'qi', 'qj', 'V', 'W']

    @accepts_compatible_units(None, None, None,
            units.elementary_charge, units.elementary_charge,
            units.kilojoules_per_mole * units.nanometers**(6),
            units.kilojoules_per_mole * units.nanometers**(12))
    def __init__(self, atom1, atom2, fudgeQQ, qi, qj, V, W):
        """
        """
        AbstractPair.__init__(self, atom1, atom2)
        self.fudgeQQ = fudgeQQ
        self.qi = qi
        self.qj = qj
        self.V = V
        self.W = W

    def get_parameters(self):
        return (self.atom1, self.atom2, self.fudgeQQ, self.qi, self.qj, self.V, self.W)

    def __repr__(self):
        return "{0} {1}: {2} {3} {4} {5}".format(self.atom1, self.atom2, self.atom3,
                self.fudgeQQ, self.qi, self.qj, self.V, self.W)

    def __str__(self):
        return "{0} {1}: {2} {3} {4} {5}".format(self.atom1, self.atom2, self.atom3,
                self.fudgeQQ, self.qi, self.qj, self.V, self.W)
