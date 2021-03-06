"""
... module:: Atom
    :platform: Unix

.. moduleuthor:: Christoph Klein <ctk3b@virginia.edu>, Christopher
Lee <ctl4f@virginia.edu>
"""
import intermol.unit as units
from converter import convert_units


class Atom(object):
    __slots__ = ['index', 'name', 'residue_index', 'residue_name',
            '_position', '_velocity', '_force', '_atomtype', 'bondtype',
            'atomic_number',
            'cgnr', '_mass', '_charge', 'ptype', '_sigma', '_epsilon']
    def __init__(self, index, name=None, residue_index=-1, residue_name=None):
        """Create an Atom object

        Args:
            index (int): index of atom in the molecule
            name (str): name of the atom (eg., N, C, H, O)
            residue_index (int): index of residue in the molecule
            residue_name (str): name of the residue (eg., THR, CYS)
        """
        self.index = index
        self.name = name
        self.residue_index = residue_index
        self.residue_name = residue_name
        self._position = [0 * units.nanometers,
                            0 * units.nanometers,
                            0 * units.nanometers]
        self._velocity = [0 * units.nanometers / units.picosecond,
                            0 * units.nanometers / units.picosecond,
                            0 * units.nanometers / units.picosecond]
        self._force = [0 * units.kilojoules_per_mole / units.nanometers,
                        0 * units.kilojoules_per_mole / units.nanometers,
                        0 * units.kilojoules_per_mole / units.nanometers]

        # These are added after data is read in and come from [ atomtypes ]
        self._atomtype = dict()
        self.bondtype = None
        self.atomic_number = None
        self.cgnr = None
        self._mass = dict()
        self._charge = dict()
        self.ptype = "A"
        self._sigma = dict()
        self._epsilon = dict()

    def getAtomType(self, index=None):
        """Gets the atomtype

        Args:
            index (str): the value corresponding with type precedence (A Type, B Type)

        Returns:
            atomtype (list, str): Returns the atomtype list or the value at
                                  index if index is specified
        """
        if index:
            return self._atomtype[index]
        return self._atomtype

    def setAtomType(self, index, atomtype):
        """Sets the atomtype

        Args:
            atomtype (str): the atomtype of the atom
            index (str): the value corresponding with type precedence (A Type, B Type)
        """
        self._atomtype[index] = atomtype



    def setSigma(self, index, sigma):
        """Sets the sigma

        Args:
            sigma (float): sigma of the atom
            index (int): index to insert at
        """
        self._sigma[index] = sigma

    def getSigma(self, index=None):
        """
        """
        if index:
            return self._sigma[index]
        return self._sigma

    def setEpsilon(self, index, epsilon):
        """Sets the epsilon

        Args:
            epsilon (float): epsilon of the atom
            index(int): index corresponding to epsilon
        """
        self._epsilon[index] = epsilon

    def getEpsilon(self, index=None):
        """
        """
        if index:
            return self._epsilon[index]
        return self._epsilon

    def setCgnr(self, index, cgnr):
        """Sets the Cgnr

        Args:
            cgnr (int): The charge group number
            index (int): the value corresponding with cgnr precedence

        """
        self._cgnr[index] = cgnr

    def getCgnr(self, index=None):
        """Gets the Cgnr

        Args:
            index (int): the index to retrieve, defaults to None

        Returns:
            cngr (dict, int): returns the index or the dictionary depending
                              on if index is set
        """
        if index:
            return self._cgnr[index]
        return self._cgnr


    def setPosition(self, x, y, z):
        """Sets the position of the atom

        Args:
            x (float): x position
            y (float): y position
            z (float): z position
        """
        unit = units.nanometers
        x = convert_units(x, unit)
        y = convert_units(y, unit)
        z = convert_units(z, unit)
        self._position = [x, y, z]

    def getPosition(self):
        """Gets the position for the atom

        Returns:
            Tuple [x, y, z]
        """
        return self._position

    def setVelocity(self, vx, vy, vz):
        """Sets the velocity of the atom

        Args:
            vx (float): x velocity
            vy (float): y velocity
            vz (float): z velocity
        """
        unit = units.nanometers / units.picoseconds
        vx = convert_units(vx, unit)
        vy = convert_units(vy, unit)
        vz = convert_units(vz, unit)
        self._velocity = [vx, vy, vz]

    def getVelocity(self):
        """Gets the velocity of the atom

        Returns:
            Tuple [vx, vy, vz]
        """
        return self._velocity

    def setForce(self, fx, fy, fz):
        """Sets the force of the atom

        Args:
            fx (float): x force
            fy (float): y force
            fz (float): z force
        """
        unit = units.kilojoules_per_mole * units.nanometers**(-1)
        fx = convert_units(fx, unit)
        fy = convert_units(fy, unit)
        fz = convert_units(fz, unit)
        self._force = [fx, fy, fz]

    def getForce(self):
        """Gets the force of the atom

        Returns:
            Tuple [fx, fy, fz]
        """
        return self._force

    def setMass(self, index, mass):
        """Sets the mass of the atom

        Args:
            mass (float): mass of the atom
            index (str): the index corresponding with mass precedence (A Mass, B Mass)
        """
        unit = units.amu
        self._mass[index] = convert_units(mass, unit)

    def getMass(self, index=None):
        """Gets the mass of the atom

        Returns:
            mass (float): mass of the atom
            index (str): index to retrieve
        """
        if index:
            return self._mass[index]
        return self._mass

    def setCharge(self, index, charge):
        """Sets the charge of the atom

        Args:
            charge (float): Charge of the atom
            index (int): the index corresponding with charge precedence
        """
        unit = units.elementary_charge
        self._charge[index] = convert_units(charge, unit)

    def getCharge(self, index=None):
        """Gets the charge of the atom

        Args:
            index (int): index of the charge to retrieve defaults to None

        Returns:
            charge (float): Charge of the atom
        """
        if index:
            return self._charge[index]
        return self._charge

    def __repr__(self):
        return 'Atom({0}, {1})'.format(self.index, self.name)

    def __cmp__(self, other):
        return self.index - other.index

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(self.index)
