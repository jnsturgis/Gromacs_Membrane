""" This module is designed to analyse membrane bilayer structures """

import md_struct

class Bilayer:
    """
    A class to contain a bilayer structure, and provides routines for analysing it.
    """
    def __init__( self, upper_leaf: MDStruct, lower_leaf: MDStruct ):
        self.upper_hg = upper_leaf
        self.lower_hg = lower_leaf
        self.energy   = 0.0
        pass

    def report( self, verbosity: int ) -> None:
        """Print a report about the membrane structure that is more or less verbose"""
        print( "***This is a membrane report.***" )
        print( 'Membrane width from difference in z position of leaflets is:' +
               f'{self.upper_hg.centroid()[2] - self.lower_hg.centroid()[2]:.3f} nm.' )

        # This next section analyses the membrane as an elastic bilayer for
        # curvature and bending constant
        # thickness and spring constant
        if verbosity > 0 :
            print( f'Upper leaflet area is {t_area:9.2f} nm$^2$, ' +
                   f'giving {t_area / float(len(p_top)):8.4f} nm$^2$/PO$_4$.' )
            print( f'Lower leaflet area is {b_area:9.2f} nm$^2$, ' +
                   f'giving {b_area / float(len(p_bottom)):8.4f} nm$^2$/PO$_4$.' )

        # Energy density
        if self.energy != 0.0 :
            print ( f'Membrane energy density is {energies[i]/(t_area+b_area):9.2f} kJ/nm².' )

        print( "" )

    def thickness( self ):
        """Analyse the membrane thickness"""
        pass
