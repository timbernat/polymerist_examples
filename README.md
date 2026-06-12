# polymerist_examples
A collection of demos, vignettes, and tutorials for the `polymerist` package for polymer structure generation and molecular dynamics (MD) simulation preparation (https://github.com/timbernat/polymerist)

## Installation
First, create a local clone of this repository in the directory of choice by running:
```sh
git clone https://github.com/timbernat/polymerist_examples
cd polymerist_examples
```
Next, follow the installation instructions in the [polymerist install docs](https://polymerist.readthedocs.io/en/latest/installation/base.html) to set up a conda environment which includes the `polymerist` package
At this point, you should be good to go and try out the examples! 

## Index of Examples
### Series 1 - [Building polymers](1-polymerization/1.0-index.ipynb)
---
* 1.1 - [Introduction: nylon polymers](1-polymerization/1.1-nylon_basics.ipynb)
* 1.2 - [Simple vinyl polymers and autopolymerization](1-polymerization/1.2-vinyl_autopolymerization.ipynb)
* 1.3 - [Kapton polyimides and with multiple intermonomer bonds](1-polymerization/1.3-polyimide_multibond_cycles.ipynb)
* 1.4 - [Crosslinkable MPD-TMC polyamides](1-polymerization/1.4-MPD-TMC_polyamides.ipynb)
* 1.5 - [PEG-PLGA block copolymers ](1-polymerization/1.5-PEG-PLGA_copolymers.ipynb)
* 1.6 - [Conjugated thiophenyl polymers with arbitrary sidechains](1-polymerization/1.6-functionalized_polythiophenes.ipynb)

### Series 2 - [Preparing systems containing polymers](2-preparation/2.0-index.ipynb)
---
* 2.1 - [Loading a polymer structure into OpenFF](2-preparation/2.1-loading_polymer_topology.ipynb)
* 2.2 - [Polymer metadata and partial charge assignment](2-preparation/2.2-preparing_individual_polymers.ipynb)
* 2.3 - [Solvation and packing of polymer melts](2-preparation/2.3-melt_packing_and_solvation.ipynb)
* 2.4 - [Reduction Charge Transfer (RCT) for generating custom library charges](2-preparation/2.4-RCT_demo.ipynb)

### Series 3 - [Running polymer simulations](3-workflows/3.0-index.ipynb)
---
* 3.1 - [Exporting polymer systems to common MD engines](3-workflows/3.1-MD_export_with_Interchange.ipynb)
* 3.2 - [Reproducibly serializing OpenMM simulations](3-workflows/3.2-serializable_simulation_parameters.ipynb)
* 3.3 - [Running series of OpenMM simulations](3-workflows/3.3-running_openmm_simulations.ipynb)
* 3.4 - [Start-to-finish polymer simulation workflow for an ATRP polymerization](3-workflows/3.4-full_workflow_demo.ipynb)

### Series 4 - [Miscellaneous](4-miscellaneous/4.0-index.ipynb)
---
* 4.1 - [Exporting arbitrary Python objects to JSON](4-miscellaneous/jsonification.ipynb)
* 4.2 - [Automated ring piercing detection (PINPRICS)](4-miscellaneous/ring_piercing.ipynb)
---

Feel free to email timotej.bernat@colorado.edu for questions or requests for more examples.