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
        atoms - with type, ID, 3D position and velocity, organized in residues with type and ID;
        box - describing a bounding box dimensions
        system_name - a description more or less intelligable of what it is.
    """
    def __init__(self):
        """
        Constructor - make a minimal valid structure (nothing in a point).
        """
        self.title = "None"
        self.n_atoms = 0
        self.atoms = []
        self.box = [0.0, 0.0, 0.0]
        pass

    def read_gro( self, filename: str ) -> None:
        """
        Read the structure that is included in the file "filename". 
        6/11/2024 Version 1 taken from above.
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
        # TODO: check here that the result is valid.
        # pass

    def write_gro( self, filename: str ) -> None:
        """
        Write a valid .gro file into the named file using the structure. Note that this
        can cause loss of precision as floats are truncated etc.
        """
        # pass

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
        # pass

    def select( self, expression: str ) -> MDStruct :
        """
        Create a new structure from a selection of the structure based on the expression.
        Currently this is performed by several specialized versions below.
        """
        new_structure = MDStruct()
        new_structure.title = "Select using "+ expression+" of "+self.title
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
        print( f"The structure contains {self.n_atoms} atoms in {self.n_residues()} residues." )
        print()
        # pass

    def match( self, expression: str, atom: list ) -> bool :
        """Does the atom match the selection expression?"""
        return False
        # pass
