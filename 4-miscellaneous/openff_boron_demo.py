'''Minimal example of how to obtain a fully-specified boron-containing molecule using OpenFF and polymerist tools'''

from pathlib import Path
import numpy as np

from openmm.unit import nanometer
from openff.toolkit import Molecule, ForceField
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper

from polymerist.mdtools.openfftools import boxvectors
from polymerist.mdtools.openfftools.topology import topology_to_sdf
from polymerist.mdtools.openmmtools.serialization import serialize_openmm_pdb
from polymerist.mdtools.openfftools.partialcharge import molchargers


# CONFIGURE PARAMETERS HERE
molname = '4-vinylphenyl-boronic_acid'
box_padding = 1*nanometer # how far beyond the tight bounding box of the polymer to extend the periodic box

# define molecule name and output directory
outdir = Path(molname)
outdir.mkdir(exist_ok=True)

# initialize molecule from SMILES string
smi = 'C=C-c1ccccc1-B(O)(O)'
offmol = Molecule.from_smiles(smi)

# assign partial charges
charger = molchargers.EspalomaCharger()
charged_mol = charger.charge_molecule(offmol)

# assign coordinates and minimal periodic box vectors
charged_mol.generate_conformers(n_conformers=1, toolkit_registry=RDKitToolkitWrapper())
charged_top = charged_mol.to_topology()

box_dims = boxvectors.get_topology_bbox(charged_top)
box_vectors = boxvectors.box_vectors_flexible(box_dims)
box_vectors = boxvectors.pad_box_vectors_uniform(box_vectors, box_padding)

# assign force field parameters
ffname = 'smirnoff99Frosst-1.0.0' # this was the only OpenFF forcefield I could find which could parameterize Boron 
ff = ForceField(f'{ffname}.offxml')
inc = ff.create_interchange(charged_mol.to_topology(), charge_from_molecules=[charged_mol])
inc.box = box_vectors

# output PDB and (the preferable) SDF Topology files and GMX files
topology_to_sdf(     outdir / f'{molname}.sdf', charged_top)
serialize_openmm_pdb(outdir / f'{molname}.pdb', inc.to_openmm_topology(), inc.positions.to_openmm())
inc.to_gromacs(prefix=f'{molname}/{molname}')