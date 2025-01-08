'''Example of how to embed coordinates for and save images of repeat unit fragments for the MPD-TMC polyamide reactions'''

import random
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdMolTransforms, AllChem

from polymerist.genutils.containers import RecursiveDict
from polymerist.genutils.sequences.seqops import int_complement


from polymerist.polymers.monomers import specification
from polymerist.rdutils.bonding import portlib, combined_rdmol
from polymerist.rdutils import rdkdraw

rdkdraw.disable_substruct_highlights()
rdkdraw.set_rdkdraw_size(400, 3/2)

from polymerist.rdutils.labeling.molwise import clear_atom_isotopes
from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from polymerist.rdutils.reactions.reactors import PolymerizationReactor


# Initialize receiving directories
RXN_DIR = Path('rxn_demo_files')
RXN_DIR.mkdir(exist_ok=True)

example_dir = RXN_DIR / 'MPD-TMC'
example_dir.mkdir(exist_ok=True)

image_dir = example_dir / 'rxn_fragment_images'
image_dir.mkdir(exist_ok=True)

# Define reactants and reaction SMARTS
reactant_smiles = {
    'MPD' : 'c1ccc(N)cc1N',                         # m-phenylenediamine
    'TMC' : 'c1c(C(=O)Cl)cc(C(=O)Cl)cc1(C(=O)Cl)',  # trimesoyl chloride
}
rxnsmarts = '[#7:1](-[*:2])(-[#1])-[#1:3].[#17]-[#6:4](=[#8:5])-[*:6]>>[#7:1](-[*:2])(-[#1:3])-[#6:4](=[#8:5])-[*:6]'

# initialize reactants
DIHEDRAL_QUERY = Chem.MolFromSmarts('*~*~*~*')

reactants : list[Chem.Mol] = []
for molname, smi in reactant_smiles.items():
    exp_smi = specification.expanded_SMILES(smi, assign_map_nums=False)
    monomer = Chem.MolFromSmiles(exp_smi, sanitize=False)
    Chem.SanitizeMol(monomer, sanitizeOps=specification.SANITIZE_AS_KEKULE)
    reactants.append(monomer)
    Draw.MolToImageFile(monomer, image_dir / f'{molname}.png')

    # conformer embedding
    errcode = AllChem.EmbedMolecule(monomer, randomSeed=0xf00d, useRandomCoords=True, ETversion=2)
    assert(errcode == 0)
    AllChem.MMFFOptimizeMolecule(monomer)
    conf = monomer.GetConformer(0)

    rdMolTransforms.CanonicalizeConformer(conf)
    Draw.MolToImageFile(monomer, image_dir / f'{molname}_ETKDG.png')

    # rotate all dihedrals randomly
    for dihed_ids in monomer.GetSubstructMatches(DIHEDRAL_QUERY):
        try:
            rdMolTransforms.SetDihedralDeg(conf, *dihed_ids, random.uniform(0, 360))
        except ValueError:
            pass
        
    rdMolTransforms.CanonicalizeConformer(conf)
    Draw.MolToImageFile(monomer, image_dir / f'{molname}_rand_torsion.png')
        
    # save to PDB file - TODO: figure out how to inject monomer names without disrupting RDKit defaults for all other AtomPDBResidueInfo
    Chem.MolToPDBFile(monomer, example_dir / f'{molname}.pdb', flavor=32)
    monomer.RemoveAllConformers()

# initialize reaction
rxn = AnnotatedReaction.from_smarts(rxnsmarts)
reactor = PolymerizationReactor(rxn_schema=rxn)

# perform fragment generation via polymerization
PORT_ID_STREAM = int_complement([0])
num_linkers_assigned : int = 0

all_dimers = {}
all_frags = RecursiveDict()
for i, (dimers, frags) in enumerate(reactor.propagate(reactants), start=1):
    for dimer in dimers:
        all_dimers[i] = dimer
        Draw.MolToImageFile(dimer, image_dir / f'dimer_step_{i}.png')

    for frag, name in zip(frags, reactant_smiles.keys()):
        valence = portlib.get_num_ports(frag)
        all_frags[name][valence] = frag
        Draw.MolToImageFile(frag, image_dir / f'{name}_{valence}-site.png')

# react remaining MPD site
reactants_post = [
    all_frags['MPD'][1],
    all_frags['TMC'][2],
]

for j, (dimers, frags) in enumerate(reactor.propagate(reactants_post), start=i+1): # resume progress
    for dimer in dimers:
        all_dimers[j] = dimer
        Draw.MolToImageFile(dimer, image_dir / f'dimer_step_{j}.png')

    for frag, name in zip(frags, reactant_smiles.keys()):
        valence = portlib.get_num_ports(frag)
        all_frags[name][valence] = frag
        Draw.MolToImageFile(frag, image_dir / f'{name}_{valence}-site.png')

# adjust port flavors
max_frags = []
for name, frag_dict in all_frags.items():
    max_frags.append(frag_dict[max(frag_dict.keys())])

combo = combined_rdmol(*max_frags, assign_map_nums=False)
for i, linker_id in enumerate(portlib.get_linker_ids(combo), start=1):
    combo.GetAtomWithIdx(linker_id).SetIsotope(i)

final_frags = Chem.GetMolFrags(combo, asMols=True)
for frag, name in zip(final_frags, reactant_smiles.keys()):
    Draw.MolToImageFile(frag, image_dir / f'{name}_final_fragment.png')