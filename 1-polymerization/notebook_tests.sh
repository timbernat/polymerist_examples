#!/bin/bash

py.test --nbval-lax --current-env -vv \
    1.1-nylon_basics.ipynb \
    1.2-vinyl_autopolymerization.ipynb \
    1.3-polyimide_multibond_cycles.ipynb \
    1.4-MPD-TMC_polyamides.ipynb \
    1.5-PEG-PLGA_copolymers.ipynb \
    1.6-functionalized_polythiophenes.ipynb \