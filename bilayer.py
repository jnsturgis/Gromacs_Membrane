""" This module is designed to analyse membrane bilayer structures """

from md_struct import MDStruct

class Bilayer:
    """
    A class to contain a bilayer structure, and provides routines for analysing it.
    """
    def __init__( self, upper_leaf: MDStruct, lower_leaf: MDStruct ):
        self.upper_hg = upper_leaf
        self.lower_hg = lower_leaf
        self.energy   = 0.0
        self.t_area   = 0.0
        self.b_area   = 0.0

    def n_headgroup( self, leaflet: int ) -> int:
        if leaflet == 0:
            return self.lower_leaf.n_residues()
        if leaflet == 1:
            return self.upper_leaf.n_residues()
        else :
            raise( ValueError, "Unknown leaflet." )

    def report( self, verbosity: int ) -> None:
        """Print a report about the membrane structure that is more or less verbose"""
        print( "***This is a membrane report.***" )
        print( 'Membrane width from difference in z position of leaflets is:' +
               f'{self.upper_hg.centroid()[2] - self.lower_hg.centroid()[2]:.3f} nm.' )

        # This next section analyses the membrane as an elastic bilayer for
        # curvature and bending constant
        # thickness and spring constant
        if verbosity > 0 :
            print( f'Upper leaflet area is {self.t_area:9.2f} nm$^2$, ' +
                   f'giving {self.t_area / float(self.n_headgroup(1))):8.4f} nm$^2$/PO$_4$.' )
            print( f'Lower leaflet area is {self.b_area:9.2f} nm$^2$, ' +
                   f'giving {self.b_area / float(self.n_headgroup(0))):8.4f} nm$^2$/PO$_4$.' )

        # Energy density
        if self.energy != 0.0 :
            print ( f'Membrane energy density is {self.energy/(self.t_area+self.b_area):9.2f} kJ/nm².' )

        print( "" )

    def thickness( self ):
        """Analyse the membrane thickness"""

