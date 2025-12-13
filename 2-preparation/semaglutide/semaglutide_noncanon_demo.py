import json
from pathlib import Path
from openff.toolkit import Topology


# define and save custom substructures for non-canonical amino acid residues
custom_substructs = {
    '2-aminoisobutyric_acid': [
        '[*:1]-[#6D3+0:2](=[#8D1+0:3])-[#6D4+0:4](-[#6D4+0:5](-[#1D1+0:9])(-[#1D1+0:10])-[#1D1+0:11])(-[#6D4+0:6](-[#1D1+0:12])(-[#1D1+0:13])-[#1D1+0:14])-[#7D3+0:7](-[*:8])-[#1D1+0:15]',
    ],
    'Amine-PEG-carboxylic_acid': [ # define middle monomer and amine-cation end group. 
        '[*:1]-[#6D3+0:2](=[#8D1+0:3])-[#6D4+0:4](-[#8D2+0:5]-[#6D4+0:6](-[#6D4+0:7](-[#8D2+0:8]-[#6D4+0:9](-[#6D4+0:10](-[#7D4+1:11](-[#1D1+0:22])(-[#1D1+0:23])-[#1D1+0:24])(-[#1D1+0:20])-[#1D1+0:21])(-[#1D1+0:18])-[#1D1+0:19])(-[#1D1+0:16])-[#1D1+0:17])(-[#1D1+0:14])-[#1D1+0:15])(-[#1D1+0:12])-[#1D1+0:13]',
        '[*:1]-[#6D3+0:2](=[#8D1+0:3])-[#6D4+0:4](-[#8D2+0:5]-[#6D4+0:6](-[#6D4+0:7](-[#8D2+0:8]-[#6D4+0:9](-[#6D4+0:10](-[#7D3+0:11](-[*:12])-[#1D1+0:23])(-[#1D1+0:21])-[#1D1+0:22])(-[#1D1+0:19])-[#1D1+0:20])(-[#1D1+0:17])-[#1D1+0:18])(-[#1D1+0:15])-[#1D1+0:16])(-[#1D1+0:13])-[#1D1+0:14]',
    ],
    'Lysine_sidechain' : [ # Lysine nitrogen "N" has one more inter-monomer site than usual. 17-carbon chain not present in the provided PDB, so is not included here
        '[*:1]-[#6D3+0:2](=[#8D1+0:3])-[#6D4+0:4](-[#7D3+0:5](-[*:6])-[#1D1+0:14])(-[#6D4+0:7](-[#6D4+0:8](-[#6D4+0:9](-[#6D4+0:10](-[#7D4+1:11](-[*:12])(-[#1D1+0:23])-[#1D1+0:24])(-[#1D1+0:21])-[#1D1+0:22])(-[#1D1+0:19])-[#1D1+0:20])(-[#1D1+0:17])-[#1D1+0:18])(-[#1D1+0:15])-[#1D1+0:16])-[#1D1+0:13]',
    ],
    'Arginine_aldehyde' : [ # Residue 36 is missing oxygen "OXT"; this is probably an (extremely hard-to-catch) error :)
        '[*:1]-[#7D3+0:2](-[#6D4+0:3](-[#6D4+0:4](-[#6D4+0:5](-[#6D4+0:6](-[#7D3+0:7](-[#6D3+0:8](-[#7D3+0:9](-[#1D1+0:22])-[#1D1+0:23])=[#7D3+1:10](-[#1D1+0:24])-[#1D1+0:25])-[#1D1+0:21])(-[#1D1+0:19])-[#1D1+0:20])(-[#1D1+0:17])-[#1D1+0:18])(-[#1D1+0:15])-[#1D1+0:16])(-[#6D3+0:11](=[#8D1+0:12])-[#1D1+0:26])-[#1D1+0:14])-[#1D1+0:13]',
    ],
}
substruct_path = Path('semaglutide_substructs.json')
with substruct_path.open('w') as file:
    json.dump(custom_substructs, file, indent=4)

# load into OpenFF using default AA residues + custom substructures
pdb_path = Path('structure.pdb')
offtop = Topology.from_pdb(pdb_path, _custom_substructures=custom_substructs)

# save to the more complete SDF format to avoid the need for filling in chemical info on each load
sdf_path = Path('semaglutide.sdf')
with sdf_path.open('w') as sdf_file:
    for offmol in offtop.molecules:
        offmol.to_file(sdf_file, file_format=sdf_path.suffix[1:])