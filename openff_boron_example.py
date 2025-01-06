'''Minimal example of how to obtain a fully-specified boron-containing molecule using OpenFF and polymerist tools'''

from pathlib import Path
import numpy as np

from openff.toolkit import Molecule, ForceField
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper

from polymerist.mdtools.openfftools.topology import topology_to_sdf
from polymerist.mdtools.openmmtools.serialization import serialize_openmm_pdb
from polymerist.mdtools.openfftools.partialcharge import molchargers


# define molecule name and output directory
molname = '4-vinylphenyl-boronic_acid'
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
conf = charged_mol.conformers[0]
conf_vals, conf_units = conf.magnitude, conf.units
box_vectors = np.diag(np.ptp(conf_vals, axis=0)) * conf_units

# assign force field parameters
ffname = 'smirnoff99Frosst-1.0.0' # this was the only OpenFF forcefield I could find which could parameterize Boron 
ff = ForceField(f'{ffname}.offxml')
inc = ff.create_interchange(charged_mol.to_topology(), charge_from_molecules=[charged_mol])
inc.box = np.diag(np.ptp(conf_vals, axis=0)) * conf_units # assigning box vectors (required for GMX >= 2020)

# output PDB and (the preferable) SDF Topology files and GMX files
topology_to_sdf(     outdir / f'{molname}.sdf', charged_mol.to_topology())
serialize_openmm_pdb(outdir / f'{molname}.pdb', inc.to_openmm_topology(), inc.positions.to_openmm())
inc.to_gromacs(prefix=f'{molname}/{molname}')