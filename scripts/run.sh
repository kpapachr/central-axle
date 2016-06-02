#!/bin/bash

PDB='../trajectories/average-4.pdb'
DCD='../trajectories/alignedtrj-ps123-ca.dcd'
REFERENCE='../trajectories/average-4.pdb'
python run_FlexE.py --pdb ${PDB} --dcd ${DCD} --reference ${REFERENCE} 
