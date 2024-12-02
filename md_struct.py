""" This module implements a structure from a molecular dynamics simulation."""
from __future__ import annotations

#
# Cell for programme development using classes
#
#
# Developing a structure class for aiding python based analysis
# Inspired by gropy written by Caizkun (?2019)
#

import os
import numpy as np
import numpy.typing as npt

TOLERANCE = 0.001				# Tolerance for 2 numbers being equal

class Atom():
    """
    Class to contain an atom and do operations on atoms. Each atom has:
        id   - a string saying what sort of particle it is.
        num  - an identification number (unique within a structure)
        pos  - the position as a numpy array of 3 floats in cartesian space
        vel  - the velocity as a numpy array of 3 floats in cartesian space
    """
    atid   : str
    num    : int
    resid  : str
    resnum : int
    pos    : npt.NDArray
    vel    : npt.NDArray

    def __init__(self, a_info: tuple(str,int), r_info: tuple(str,int),
                       a_pos : list[float], a_vel : list[float] ):
        """Constructor, use passed values to make an atom."""
        self.atid   = a_info[0]
        self.num    = a_info[1]
        self.resid  = r_info[0]
        self.resnum = r_info[1]
        assert len(a_pos) == 3
        assert len(a_vel) == 3
        self.pos    = np.array(a_pos)
        self.vel    = np.array(a_vel)

    def distance(self, location: npt.NDArray) -> float:
        """Find the distance of an atom from an arbitrary location"""
        assert len(location) == 3
        vec = self.pos - location
        return np.sqrt( np.sum ( vec * vec ))

    def speed( self ) -> float:
        """Find the speed of an atom from the velocity vector"""
        return np.sqrt( np.sum ( self.vel * self.vel ))

    def __eq__( self, other: Atom ) -> bool:
        """
        Are 2 atoms at the same place and same velocity and of the same type? Note they are still
        considered the same if they have different unique ID numbers and belong to residues with
        different names or numbers.
        """
        # pylint: disable=chained-comparison
        return ((self.atid == other.atid)
        	and (self.distance( other.pos ) < TOLERANCE )
        	and TOLERANCE > np.sqrt(np.sum( (self.vel - other.vel) * (self.vel - other.vel)))
        	)
class Residue():
    """
    A residue class for handling residues within structures. Each residue has:
        name  - which describes the type of residue, for example DOPC for di-oleyl phosphatidyl
                choline
        idnum - a unique id number (there can be more than one copy of a given type of residue).
                Within a valid structure no two residues can have the same unique id number.
        atoms - a list of atoms.
    """
    name  : str
    idnum : int
    atoms : list[Atom]

    def __init__(self, r_name: str, r_id: int ):
        """Create a new residue with the given name and idnumber but no atoms"""
        self.name = r_name
        self.idnum = r_id
        self.atoms = []

    def add_atom( self, atom: Atom ) -> None:
        """Add an atom to the residue"""
        self.atoms.append(atom)

    def __eq__(self, other: Residue ) -> bool:
        """Two redidues are equal if they are of the same type and their atoms are equal."""
        # pylint: disable=multiple-statements
        if self.name != other.name : return False
        if len(self.atoms) != len(other.atoms) : return False
        for atom1, atom2 in zip(self.atoms, other.atoms):
            if atom1 != atom2: return False
        return True

class MDStruct():
    """
    A structural configuration class, for example from a molecular dynamics or structure file.
    This object contains minimally:
        title    - some label describing the structure
        filename - where this information is stored or came from
        atoms    - with type, ID, 3D position and velocity, organized in residues with type and ID
        resids   - the residues in the structure with name, ID, and atoms.
        box      - describing a bounding box dimensions
    """
    title    : str
    filename : str
    atoms    : list[ Atom ]                      # Actually a list of atom structures
    resids   : list[ any ]                       # TODO Actually a list of residues
    box      : list[ float ]

    def __init__(self):
        """
        Constructor - make a minimal valid structure (nothing in a point).
        """
        self.title = "None"
        self.filename = ""
        self.atoms = []
        self.resids = []
        self.box = [0.0, 0.0, 0.0]
        pass

    def n_atoms( self ) -> int:
        """How many atoms in the structure"""
        return len(self.atoms)

    def n_resid( self ) -> int:
        """How many residues in the structure"""
        return len(self.resids)

    def is_valid( self ) -> bool:
        """Check that the structure is sane"""
        # TODO all atoms in box
        # TODO atom - residue concordance
        return  (len(self.title) > 0
               and len(self.box) == 3)

    def copy( self ) -> MDStruct:
        """Make a copy of a structure, so one can be modified and the other remains the same."""
        new = MDStruct()
        new.title = self.title                   # TODO Probably should do a deeper copy here
        new.filename = ""
        for atom in self.atoms :
            new.atoms.append( atom )             # TODO Probably should do a deeper copy here
        for resid in self.resids :
            new.resids.append( resid )           # TODO Probably should do a deeper copy here
        new.box = self.box                       # TODO Probably should do a deeper copy here
        return new

    def merge( self, other: MDStruct ) -> None:
        """
        Insert the residues and atoms of the second structure into the first, which is renamed
        in consequence.
        """
        # pylint: disable=multiple-statements
        for atom in other.atoms: self.atoms.append(atom)
        for residue in other.resids: self.resids.append(residue)
        part1 = self.filename if len(self.filename) > 0 else self.title
        part2 = other.filename if len(other.filename) > 0 else other.title
        self.title = f"{part1} plus {part2}"
        self.filename = ""

    def read_gro( self, filename: str ) -> None:
        """
        Read the structure that is included in the file "filename".
        6/11/2024 Version 1 adapted from Gropy/Gro.py on github.
        """
        # TODO files can contain multiple structures should be able to read the n'th
        # record and fail elegantly if there is not one.
        self.filename = filename
        with open(filename, 'r', encoding="utf-8") as file_id:
            for line_count, line in enumerate(file_id):
                line = line.rstrip()
                if   line_count == 0 :
                    self.title = line
                elif line_count == 1 :
                    n_atoms = int(line)
                    last_line = n_atoms + 1
                elif line_count <= last_line :
                    if len(line) > 44 :
                        at_vx = float(line[44:52])
                        at_vy = float(line[52:60])
                        at_vz = float(line[60:68])
                    else :
                        at_vx = 0.0
                        at_vy = 0.0
                        at_vz = 0.0
                    self.atoms.append( Atom( (line[5:10].strip(), int(line[0:5])),
                               (line[10:15].strip(), int(line[15:20])),
                               [ float(line[20:28]), float(line[28:36]), float(line[36:44]) ],
                               [ at_vx, at_vy, at_vz ]))
                else :
                    self.box = [ float( size ) for size in line.split() ]
        file_id.close()
        # TODO setup the residues
        current_residue = -1
        residue = None
        for atom in self.atoms:
            if atom.resnum > current_residue :
                current_residue = atom.resnum
                residue = Residue( atom.resid, atom.resnum)
                self.resids.append(residue)
            residue.add_atom(atom)
        # Check here that the result is valid.
        assert self.is_valid()
        pass

    def write_gro( self ) -> None:
        """
        Write a valid .gro file into the named file using the structure. Note that this
        can cause loss of precision as floats are truncated etc.
        Version 1.0 27/11/24 adapte from gropy/Gro.py
        """
        with open(self.filename, 'w', encoding="utf-8") as file_id:
            file_id.write(f"{self.title}\n")
            file_id.write(f" {self.n_atoms():d}\n")
            for atom in self.atoms:
                file_id.write(f"{atom.num:5d}{atom.atid:<5s}{atom.resid:>5s}{atom.resnum:5d}"
                        f"{atom.pos[0]:8.3f}{atom.pos[1]:8.3f}{atom.pos[2]:8.3f}"
                        f"{atom.vel[0]:8.4f}{atom.vel[1]:8.4f}{atom.vel[2]:8.4f}\n")
            file_id.write(f"{self.box[0]:10.5f}{self.box[1]:10.5f}{self.box[2]:10.5f}\n")
            file_id.close()

    def write_gro_as( self, filename : str) -> None:
        """Change the filename associated with a structure and then save."""
        self.filename = filename
        self.write_gro()

    def renumber(self) -> None:
        """
        Renumber atoms and residues (change the atom and residue id's) so that they are sequential.
        """
        last_resid  = -1
        curr_resid  = 0
        curr_atomid = 0
        for atom in self.atoms:
            curr_atomid += 1
            if atom[0] != last_resid :
                curr_resid += 1
            atom[0] = curr_resid
            atom[3] = curr_atomid

        assert curr_atomid == self.n_atoms()
        pass

    def select( self, expression: str ) -> MDStruct :
        """
        Create a new structure from a selection of the structure based on the expression.
        Currently this is performed by several specialized versions below.
        """
        new_structure = MDStruct()
        new_structure.title = "Select using "+ expression +" of "+ self.title
        new_structure.box   = self.box

        for atom in self.atoms:
            if new_structure.match( expression, atom ):
                new_structure.atoms.append( atom )
        new_structure.renumber()
        return new_structure

    def select_hg( self ) -> MDStruct :
        """
        Select headgroup atoms from a structure.
        """
        new_structure = MDStruct()
        new_structure.title = "Headgroups (PO4, PA4, PB4) extracted from "+self.title
        new_structure.box   = self.box

        for atom in self.atoms:
            at_name = atom[2]
            if at_name in ("PO4", "PA4", "PB4"):
                new_structure.atoms.append( atom )

        new_structure.renumber()
        return new_structure

    def select_id( self, start: int, end: int ) -> MDStruct :
        """
        Select atoms based on a range of atom_id numbers.
        """

        if start < 0 :
            start = self.n_atoms() + 1 + start
        if end   < 0 :
            end   = self.n_atoms() + 1 + end

        assert start <= end , "The 'start' must be less than 'end'"

        new_structure = MDStruct()
        new_structure.title = f"Select {start} to {end} of " + self.title
        new_structure.box   = self.box

        for atom in self.atoms:
            if atom[3] <= end and atom[3] >= start :
                new_structure.atoms.append( atom )

        new_structure.renumber()
        return new_structure

    def centroid( self ) :
        """Calculate the center of the atoms in a structure"""
        center = np.array([0.0, 0.0, 0.0])
        count = 0
        for atom in self.atoms:
            center += np.array( [atom[4], atom[5], atom[6]] )
            count += 1
        center /= count
        return center

    def report( self, verbosity: int ) -> None:
        """Print a report on the structure, depending on the level of verbosity desired."""
        print( f"This is a report on the structure:\n{self.title}" )
        if verbosity > 0 :
            print( f"The structure contains {self.n_atoms()} atoms in {self.n_resid()} residues." )
            print()
        pass

    def match( self, expression: str, atom: list ) -> bool :
        """Does the atom match the selection expression?"""
        if atom[0] > 1 :
            print( expression )
        return False

# Dunder comparison methods that make sense only comparing 2 structures
    def __eq__( self, other: MDStruct ) -> bool:
        """
        Check if two structures are equal that is have the same number of atoms and residues,
        and the atoms are of the same type with the same position and velocity (within rounding
        errors). This does not require the same file and title.
        TODO the residues should match as well in type, and atomic composition.
        """
        result = True
        result &= (self.n_atoms() == other.n_atoms())
        result &= (self.n_resid() == other.n_resid())
        if not result:
            return result
        for atom1, atom2 in zip(self.atoms, other.atoms):
            result &= (atom1 == atom2)
            if not result:
                break
        result &= np.max(np.abs(np.array(self.box) - np.array(other.box))) < TOLERANCE
        return result

    def __ne__( self, other: MDStruct ) -> bool:
        """Not equal is obvious"""
        return not self == other

    def __str__( self ) -> str:
        """Print a sensible description of the MDStruct object"""
        return f"MDStruct: '{self.title}' with {self.n_atoms()} atoms in a box."

# Test routines to check the class works OK.

def test_start() -> None:
    """Notify start of testing."""
    print("************************************************")
    print("***** Starting tests on the class MDStruct *****")

def test_end() -> None:
    """Notify end of testing."""
    print("***** Tests on MDStruct successfully done. *****")
    print("************************************************")

def test_read_write() -> None:
    """Test the file reading and writing routines"""
    errors = 0

    # Remove old temporary files
    # pylint: disable=multiple-statements
    try: os.remove("test_out.gro")
    except OSError: pass
    try: os.remove("test_out2.gro")
    except OSError: pass

    print("*                                              *")
    print("* Testing file reading and writing routines.   *")
    structure = MDStruct()                            # Create a structure with data from the file
    structure.read_gro( "test/data/DOPC2EC_100a_S6.gro" )
    print("* Successfully read test structure gro file.   *")
    if structure.is_valid() :
        print("* The structure read from the file is valid.   *")
    else :
        print("* FAILURE: An invalid structure has resulted,  *")
        print("* It will be written to test_out.gro to check. *")
        errors += 1
    structure.write_gro_as( "test_out.gro" )
    print("* Successfully wrote test structure gro file.  *")
    # todo: check that it is as it should be
    new_struct = MDStruct()
    new_struct.read_gro( "test_out.gro" )
    if new_struct == structure :
        print("* Write and read did not modify the structure  *")
    else :
        print("* FAILURE: The structure has been changed by a *")
        print("* write read cycle. Compare 'test_out.gro' and *")
        print("* 'test_out2.gro'.                             *")
        new_struct.write_gro_as( 'test_out2.gro' )
        errors += 1
    # If no errors found remove temporary files
    if not errors:
        print("* Removing temporary files as all well in I/O  *")
        try: os.remove("test_out.gro")
        except OSError: pass
    print("*                                              *")

# Unit tests that are run if the file is called directly.
if __name__ == "__main__":
    test_start()
    test_read_write()
    # test_select()
    # test_merge()
    # test_report()
    test_end()
