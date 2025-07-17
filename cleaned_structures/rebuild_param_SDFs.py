'''
For quickly regenerating parameterized SDF files of all loadable PDB structures with monomer fragments available
NOTE: !!NOT A DEMO!! Running this is pretty expensive! Only do this for recovery purposes 
'''

# Supressing annoying warnings (!must be done first!)
import warnings
warnings.catch_warnings(record=True)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import logging
logging.basicConfig(level=logging.INFO)

from collections import defaultdict
from pathlib import Path

from openff.toolkit import Topology
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
from openff.toolkit.utils.exceptions import UnassignedChemistryInPDBError, IncorrectNumConformersWarning
warnings.filterwarnings('ignore', category=IncorrectNumConformersWarning)

from polymerist.polymers.monomers.repr import MonomerGroup, LOGGER as monomer_logger
monomer_logger.setLevel(logging.WARNING) # suppress verbose fragment registration messages to allow other info to breathe

from polymerist.mdtools.openfftools.topology import topology_to_sdf, get_largest_offmol
from polymerist.mdtools.openfftools.partition import partition
from polymerist.mdtools.openfftools.partialcharge.molchargers import NAGLCharger, EspalomaCharger


overwrites_allowed : bool = True

# locate and initialize directories
# master_dir = Path('cleaned_structures')
master_dir = Path('.')
pdb_dir  = master_dir / 'pdbs'
mono_dir = master_dir / 'monos'

sdf_dir  = master_dir / 'SDF'
sdf_dir.mkdir(exist_ok=True)

# collate available monomer
charger = EspalomaCharger()
errors = defaultdict(list)
monos_available : dict[str, MonomerGroup] = {}

for mono_path in mono_dir.glob('**/*.json'):
    try:
        monogrp = MonomerGroup.from_file(mono_path)
        if isinstance(monogrp, dict): # missing class tags lead to JSON output being decoded as a dict...
            monogrp = MonomerGroup(**monogrp) # ...instead of the target type (silent error), so need to correct
            if overwrites_allowed:
                monogrp.to_file(mono_path) # rewrite file to be format-compliant
    except (TypeError, ValueError): # confusingly, the "Invalid SMARTS" ValueError is caught and reraised as a TypeError by the TypeSerializer Decoder
        errors['Invalid residue SMARTS'].append(mono_path.stem)

    monos_available[mono_path.stem] = monogrp

# load PDBs which have correcsponding monomer templates defined
pdb_candidates = [pdb_path for pdb_path in pdb_dir.glob('**/*.pdb')] # NOTE: unpacking as list solely to get headcount on number of PDBs
n_pdbs_total = len(pdb_candidates)
n_successful : int = 0

for i, pdb_path in enumerate(pdb_candidates, start=1):
    molname = pdb_path.stem
    logging.info(f'Examining PDB {i}/{n_pdbs_total} ({molname})...')
    if pdb_path.stem not in monos_available:
        logging.warning(f'Skipping {molname}...')
        errors['Missing residue definitions'].append(molname)
        continue # skip PDBs with no fragments defined
    else:
        monogrp = monos_available[molname]

    # Load residue info via OpenFF
    try:
        offtop = Topology.from_pdb(pdb_path, _custom_substructures=monogrp.monomers)
    except UnassignedChemistryInPDBError:
        logging.warning(f'Skipping {molname}...')
        errors['Bad residue cover'].append(molname)
        continue

    try:
        was_partitioned = partition(offtop)
        if not was_partitioned:
            raise ValueError
    except (ValueError, AssertionError):
        logging.warning(f'Skipping {molname}...')
        errors['Partition failed'].append(molname)
        continue

    # extract molecule, assign partial charges
    offmol = get_largest_offmol(offtop)
    cmol = charger.charge_molecule(offmol)

    # save to SDF
    ## set output SDF path to have same internal file structure as PDB dir, but relative to SDF dir and with approriate extension
    sdf_path = Path(sdf_dir, pdb_path.relative_to(pdb_dir)).with_suffix('.sdf') 
    sdf_path.parent.mkdir(exist_ok=True, parents=True) # initialize any missing parent directories
    topology_to_sdf(sdf_path, cmol.to_topology(), toolkit_registry=RDKitToolkitWrapper())

    # record success
    logging.info(f'Successfully generated parameterized SDF for {molname}!')
    n_successful += 1

logging.info(f'SDF output complete! {n_successful}/{n_pdbs_total} SDFs written without error; errors encountered are as follow:')
print(dict(errors))