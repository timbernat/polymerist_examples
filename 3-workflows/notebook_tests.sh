#!/bin/bash

py.test --nbval-lax --current-env -vv \
    3.1-MD_export_with_Interchange.ipynb \
    3.2-serializable_simulation_parameters.ipynb \
    3.3-running_openmm_simulations.ipynb \
    3.4-full_workflow_demo.ipynb \