# polymerist_examples
A collection of demos, vignettes, and tutorials for the polymerist polymer structure and MD package (https://github.com/timbernat/polymerist)

## Installation from GitHub with Conda
First create a local clone of this repository and, within it, a clone of the in-house toolkit this repo depends on:
the toolkit, first download and install . Then, run:
```sh
git clone https://github.com/timbernat/polymerist_examples
cd polymerist_examples
git clone https://github.com/timbernat/polymerist
```

Next, recreate the conda environment (after installing the [Anaconda Distribution](https://www.anaconda.com/download)) from the packaged requirements:
```sh
conda env create -n polymerist_env -f reqs.yml
conda activate polymerist_env
```
