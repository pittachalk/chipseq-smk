# Skeleton for the Snakemake tutorial

This repository hosts the skeleton code needed for the [Snakemake tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html).

To create the environment:
```
conda env create --name snakemake-tutorial --file environment.yaml
```

Activate the environment, and deactivating it:
```
conda activate snakemake-tutorial
conda deactivate
```

Chances are you might need to run `source .bashrc` to get `conda init` to work (at least that was the case for me).
