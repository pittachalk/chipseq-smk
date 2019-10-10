# Skeleton for the Snakemake tutorial

This repository hosts the skeleton code needed for the [Snakemake tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html).

## Setting up
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


## More advanced stuff
Executing as a dry run `-n`:
```
snakemake -np mapped_reads/{A,B}.bam
```

Show the DAG of jobs
```
snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg
```

Execute the rule `all` at the top of Snakefile:
```
snakemake
```

Force a specific step to run, you can add ``--reason`` to know why Snakemake chooses to run each step.
```
snakemake --forcerun samtools_sort
```


