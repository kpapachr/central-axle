from MDAnalysis import *
from prody import *
import numpy as np
import pandas as pd
import hamiltonian
import argparse

def main():

    # No log messages from prody
    prody.confProDy(verbosity='none')

    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate MDENM energies from a pdb \
                                    will calculate energy using modes from pdb\
                                    and then from reference--> crystal should\
                                    be the reference')
    parser.add_argument('--pdb', required=True,
                        help='Molecule we want to examine. It will also be used as topology')
    parser.add_argument('--dcd', required=True,
                        help='Molecule we want to examine. It will also be used as topology')
    parser.add_argument('--reference', required=True,
                        help='Rerence pdb to which we will rmsd everything')
    args = parser.parse_args()


    # Load the structures
    # a REFERENCE for the ENM and a DCD trajectory to be analysed
    pdb = parsePDB( args.pdb )
    calphas = pdb.select('calpha')
    
    # Get first the reference structure
    ref = prody.parsePDB( args.reference )
    ref_alpha = ref.select( 'calpha' )

    # then read the DCD file
    dcd = DCDFile( args.dcd )

    # Associate that with the reference structure
    dcd.setCoords(ref_alpha)
    #dcd.link(ref_alpha)

    # Loop over trajectory frames
    listFrames = list()
    listFlexE  = list() 
    for frameIndx, frame in enumerate( dcd ):
	    # Make sure we are in same reference set
	    #t = calcTransformation(calphas,ref_alpha)
	    #t.apply(calphas)

	    native = ref_alpha.copy()
	    pred = frame

	    h = hamiltonian.EDENMHamiltonian( native.getCoords() )
	    Forw_E_ED = h.evaluate_energy( pred.getCoords())

	    number_of_residues = ref.numResidues()

	    listFrames.append( frameIndx )
            listFlexE.append( Forw_E_ED / number_of_residues )

    # Write out a CSV file with the results
    df = pd.DataFrame( {'FrameIndx':listFrames, 'FlexE':listFlexE} ).set_index('FrameIndx')
    df.to_csv('FlexE-ps123.csv')

if __name__ == '__main__':
    main()

