#!/bin/bash

py.test --nbval-lax --current-env -vv \
    2.1-loading_polymer_topology.ipynb \
    2.2-preparing_individual_polymers.ipynb \
    2.3-melt_packing_and_solvation.ipynb \
    2.4-RCT_demo.ipynb \