""" This module implements a structure from a molecular dynamics simulation."""
from __future__ import annotations

#
# Cell for programme development using classes
#
#
# Developing a structure class for aiding python based analysis
# Inspired by gropy written by Caizkun (?2019)
#

import numpy as np

class MDStruct():
    """
    A structural configuration class, for example from a molecular dynamics or structure file.
    This object contains minimally: 
        title   - some label describing the structure
        n_atoms - the number of atoms in the structure
        n_resid - the number of residues in the structure
        atoms   - with type, ID, 3D position and velocity, organized in residues with type and ID
        box     - describing a bounding box dimensions
    """
    title   : str
    n_atoms : int
    atoms   : list[ any ]                        # Actually a list of atom structures
    n_resid : int                                # TODO
    resids  : list[ any ]                        # TODO Actually a list of residues
    box     : list[ float ]

    def __init__(self):
        """
        Constructor - make a minimal valid structure (nothing in a point).
        """
        self.title = "None"
        self.n_atoms = 0
        self.atoms = []
        self.n_resid = 0
        self.resids = []
        self.box = [0.0, 0.0, 0.0]
        pass

    def is_valid( self ) -> bool:
        """Check that the structure is sane"""
        # TODO all atoms in box
        # TODO atom - residue concordance
        return  (len(self.title) > 0
               and len(self.atoms) == self.n_atoms
               and len(self.resids) == self.n_resid
               and len(self.box) == 3)

    def read_gro( self, filename: str ) -> None:
        """
        Read the structure that is included in the file "filename". 
        6/11/2024 Version 1 adapted from Gropy/Gro.py.
        """
        with open(filename, 'r', encoding="utf-8") as file_id:
            for line_count, line in enumerate(file_id):
                line = line.rstrip()
                if   line_count == 0 :
                    self.title = line
                elif line_count == 1 :
                    self.n_atoms = int(line)
                    last_line = self.n_atoms + 1
                elif line_count <= last_line :
                    if len(line) > 44 :
                        at_vx = float(line[44:52])
                        at_vy = float(line[52:60])
                        at_vz = float(line[60:68])
                    else :
                        at_vx = 0.0
                        at_vy = 0.0
                        at_vz = 0.0
                    self.atoms.append([ int(line[0:5]), line[5:10].strip(),
                                        line[10:15].strip(), int(line[15:20]),
                                        float(line[20:28]), float(line[28:36]), float(line[36:44]),
                                        at_vx, at_vy, at_vz ])
                else :
                    self.box = [ float( size ) for size in line.split() ]
        file_id.close()
        # TODO setup the residues
        # Check here that the result is valid.
        assert self.is_valid()
        pass

    def write_gro( self, filename: str ) -> None:
        """
        Write a valid .gro file into the named file using the structure. Note that this
        can cause loss of precision as floats are truncated etc.
        Version 1.0 27/11/24 adapte from gropy/Gro.py
        """
        with open(filename, 'w', encoding="utf-8") as file_id:
            file_id.write(f"{self.title}\n")
            file_id.write(f" {self.n_atoms:d}\n")
            for atom in self.atoms:
                file_id.write(f"{atom[0]:5d}{atom[1]:<5s}{atom[2]:>5s}{atom[3]:5d}"
                        f"{atom[4]:8.3f}{atom[5]:8.3f}{atom[6]:8.3f}"
                        f"{atom[7]:8.4f}{atom[8]:8.4f}{atom[9]:8.4f}\n")
            file_id.write(f"{self.box[0]:10.5f}{self.box[1]:10.5f}{self.box[2]:10.5f}\n")
            file_id.close()
        pass

    def n_residues( self ) -> int:
        """
        Count the number of residues in the structure.
        This requires and assumes that all atoms of an individual residue are grouped together,
        it also assumes that negative residue id's are invalid.
        """
        count = 0
        last_resid = -1
        for atom in self.atoms:
            if atom[0] != last_resid :
                last_resid = atom[0]
                count += 1
        return count

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

        assert curr_atomid == self.n_atoms
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
                new_structure.n_atoms += 1
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
                new_structure.n_atoms += 1
                new_structure.atoms.append( atom )

        new_structure.renumber()
        return new_structure

    def select_id( self, start: int, end: int ) -> MDStruct :
        """
        Select atoms based on a range of atom_id numbers.
        """

        if start < 0 :
            start = self.n_atoms + 1 + start
        if end   < 0 :
            end   = self.n_atoms + 1 + end

        assert start <= end , "The 'start' must be less than 'end'"

        new_structure = MDStruct()
        new_structure.title = f"Select {start} to {end} of " + self.title
        new_structure.box   = self.box

        for atom in self.atoms:
            if atom[3] <= end and atom[3] >= start :
                new_structure.n_atoms += 1
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
            print( f"The structure contains {self.n_atoms} atoms in {self.n_residues()} residues." )
            print()
        pass

    def match( self, expression: str, atom: list ) -> bool :
        """Does the atom match the selection expression?"""
        if atom[0] > 1 :
            print( expression )
        return False

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
    print("*                                              *")
    print("* Testing file reading and writing routines.   *")
    structure = MDStruct()                            # Create a structure with data from the file
    structure.read_gro( "test/data/DOPC2EC_100a_S6.gro" )
    structure.is_valid()
    print("* Successfully read test structure gro file.   *")
    structure.write_gro( "test_out.gro" )
    print("* Successfully wrote test structure gro file.  *")
    # todo: check that it is as it should be
    # todo: write and then reread and check the same.
    print("*                                              *")

# Unit tests that are run if the file is called directly.
if __name__ == "__main__":
    test_start()
    test_read_write()
    # test_select()
    # test_merge()
    # test_report()
    test_end()
