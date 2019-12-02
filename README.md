# chipseq-snakemake
This is a custom built pipeline to perform ChIP-seq analysis, developed using the Snakemake workflow manager. There is no confidental data on this repo, so I welcome anyone to use the pipeline for their own analyses. Please report any issues or bugs. The pipeline executes the following steps:

* Concatenates FASTQ files from different lanes.
* Trimming with Trimmomatic.
* FASTQC quality control.
* Map to reference genome with BWA.
* Peak calling with MACS2.
* Find common peaks between each pair of replicates (sample peaks).
* Find the union of sample peaks in the whole dataset.
* Plot heatmaps of the ChIP signal over these peaks (by sample).
* Perform PCA of all replicates of all samples.

## Contents of this repo
The most important components are:

* `Snakefile`: rules for each step of the pipeline.
* `config.yaml`: configuration file with info about samples (directory, names).
* `envs/`: configuration files specifying package dependencies for the Conda package manager.

Less important is `cluster.json`, which contains custom instructions to utilise multiple cores on a LSF-based machine (like the one I have access to). `data/` contains sample data, and `adapter` provides adapter sequences for trimming.

## Setup
Install the Conda package manager ([link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)). My current version is 4.7.12. Also, I recommend not letting the installer overwriting your default path for Python.

Next, use Conda to build the Python environments in `envs/`. `py3.yml` specifies a Python 3 environment that works for all packages, except for MACS2 which must be executed in Python 2 (`py2.yml`).
```
# create environments
conda env create -f envs/py3.yml
conda env create -f envs/py2.yml

# verify that these have been installed
conda env list
conda info --envs  # alternatively
```

## Executing the pipeline
Activate the Python3 environment.
```
conda activate chipseq
```

Execute Snakemake as a dry run first. This should print all the commands. Make sure everything works, then rerun the command with `-p` replacing `-np`.
```
snakemake -np --use-conda
```

You can also generate a diagram showing the directed acyclic graph of the workflow.
```
snakemake --dag | dot -Tsvg > dag.svg
```

If you're on a SLURM cluster, you can use the settings in `cluster.json` to utilise multiple cores:
```
# create a directory for slurm log messages
mkdir slurm/

# run snakemake
snakemake --cores 8 --restart-times 1 \
  --printshellcmds --keep-going --use-conda \
  --cluster-config ./cluster.json --cluster \
  "sbatch -J {cluster.jobname} -o {cluster.output} \
  -c {cluster.cores} -N {cluster.nodes} -p {cluster.partition} \
  --mail-user {cluster.mailuser} --mail-type {cluster.mailtype}"
```

### Disclaimer about package specific bugs
#### MACS2 (and why this repo does not work out of the box)
There is a weird bug in the peak caller MACS2 which causes it to crash if the sample BAM files contains no peaks in it. This is documented [here](https://github.com/taoliu/MACS/issues/108), and is beyond my capability to fix, especially in the context of a Conda environment.

For most intents and purposes, this error rarely comes up in real data. If the pipeline crashes at the MACS2 step and throws this error, you should interpret this as the absence of peaks (and move along...).

This is indeed the case for the sample data in this repo, but running things as a dry run `-n` should still work.

#### IDR
IDR does not work if there are fewer than 20 common peaks between replicates.

## Optional reading (mostly for my records)
### Wishlist and to-dos
* Add differential binding analysis, based on [here](https://github.com/taoliu/MACS/wiki/Call-differential-binding-events)
* Simulate real data for chip seq peak calling (so this repo works out of the box)

### More advanced options for Snakemake 
Link to the original Snakemake [tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html), and how to [manage environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

Snakemake operates on a file basis, meaning that you can run something like `snakemake <outputfilename>`, and it tries to identify the rules needed to generate the file(s). Otherwise, it attempts to generate the files specified by the rule `all` (first rule in `Snakefile`).

Snakemake only reruns steps in a pipeline if necessary. Use `--force` or `--forceall` to alter that behaviour. You can add ``--reason`` to know why Snakemake chooses to run each step.

`snakemake --forcerun <rule>` forces a rule to be run.

### Environments
In case you change versions of packages, you can update the `yml` files by exporting your environment:
```
conda env export > environment.yml
```

Or create an environment from scratch (after which you install packages using `conda install`):
```
conda env create --name <envname> --file environment.yaml
```

To remove a Conda environment:
```
conda remove --name <envname> --all
```

Chances are you might need to run `source .bashrc` to get `conda init` to work (at least that was the case for me).

### Quirks in Snakemake 
* If the `conda` directive in a rule has spaces instead of tabs, it does not work.
* Directives like `params` and `threads` are fixed. Custom directives are not allowed.
* Config files can be interpreted as nested lists. `config["hey"]["oh"]`

#### expand
I recommend using `expand` only in the rule `all`. According to the docs, `expand` is run during the initialisation phase, during which wildcard values are not known. Therefore, the `"{outdir}{sample}.out"` values must be taken from the config file.

I can't seem to get `expand` to be robust in other rules. Bug messages always revolve around the same command being called twice, which might happen due to the nature of our analyis requiring the concatenation of two files.  

My workaround is to extract the config parameters at the start of the file, rather than using expand later on. This can be accessed with the simple `+` function (i.e. `outdir + "{sample}_joined.txt"`). Seems clunky, but it works for me.

### Minimal skeleton of Snakemake file for reference
```
configfile: "config.yaml"
indir  = config["indir"]
outdir = config["outdir"]

rule all:
	input:
	    expand("{outdir}{sample}.out", sample=config["samples"], outdir = config["outdir"])

rule runanalysis:
	input:
		lambda x: map(lambda y: indir + y, config["samples"][x.sample])
	output:
		outdir + "{sample}_joined.txt"
	shell:
		"cat {input} > {output}"
```

And it's corresponding `config.yaml` file:
```
samples:
  A: 
    - A_1.txt
    - A_2.txt
  B:
    - B_1.txt
    - B_2.txt
    - B_3.txt

indir: input/

outdir: output/
```

Having `print(config)` in the Snakefile returns a dictionary of variables:
```
{'samples': {'A': ['A_1.txt', 'A_2.txt'], 'B': ['B_1.txt', 'B_2.txt', 'B_3.txt'],}, 'indir': 'input/', 'outdir': 'output/'}
```

##### {sample}
`{sample}` that outside the rule `all` are unusual. Intuitively, one might expect is that it's the keys in `config["samples]` dictionary, but it is NOT. Rather, it is a **wildcard** that Snakemake attempts to match from the output file specified. It seems to be equivalent to `wildcards.sample` when used in a shell directive.

##### a list of lists
`lambda x: map(lambda y: indir + y, config["samples"][x.sample])` is ugly, but what it's 

`config["samples] = {'A': ['A_1.txt', 'A_2.txt'], 'B': ['B_1.txt', 'B_2.txt', 'B_3.txt'],}`

* Snakemake uses the outer `lambda x:` to apply a function to a list of lists: `[['A_1.txt', 'A_2.txt'], ['B_1.txt', 'B_2.txt', 'B_3.txt']]`.
* That function is to append the indir to each sample, without opening up the list: `[['input/A_1.txt', 'input/A_2.txt'], ['input/B_1.txt', 'input/B_2.txt', 'input/B_3.txt']]`. We use `map` for this.
* When the shell command is run `"cat {input} > {output}"`, the elements of the list are joined with a space. `cat input/A_1.txt input/A_2.txt > output/A_joined.txt`

