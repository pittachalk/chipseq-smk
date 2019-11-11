# Skeleton for the Snakemake tutorial

## Quirks I noticed about Snakemake
### Why are my variables not recognised??!!
Here is a (relatively) minimal example of a Snakemake file, representative of my workflow/logic.
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

#### expand
I recommend using `expand` only in the rule `all`. According to the docs, `expand` is run during the initialisation phase, during which wildcard values are not known. Therefore, the `"{outdir}{sample}.out"` values must be taken from the config file.

I can't seem to get `expand` to be robust in other rules. Bug messages always revolve around the same command being called twice, which might happen due to the nature of our analyis requiring the concatenation of two files.  

My workaround is to extract the config parameters at the start of the file, rather than using expand later on. This can be accessed with the simple `+` function (i.e. `outdir + "{sample}_joined.txt"`). Seems clunky, but it works for me.

#### {sample}
`{sample}` that is not in the rule `all` are weird. I do not know where it's getting the values from. My intuition is that it's the keys in `config["samples]` dictionary, but I can't confirm this. Note that this `{sample}` is not equivalent to that in the `expand` command in the rule `all`. It seems to be equivalent to `wildcards.sample` when used in a shell.

#### a list of lists
`lambda x: map(lambda y: indir + y, config["samples"][x.sample])` is ugly, but what it's 

`config["samples] = {'A': ['A_1.txt', 'A_2.txt'], 'B': ['B_1.txt', 'B_2.txt', 'B_3.txt'],}`

* Snakemake uses the outer `lambda x:` to apply a function to a list of lists: `[['A_1.txt', 'A_2.txt'], ['B_1.txt', 'B_2.txt', 'B_3.txt']]`.
* That function is to append the indir to each sample, without opening up the list: `[['input/A_1.txt', 'input/A_2.txt'], ['input/B_1.txt', 'input/B_2.txt', 'input/B_3.txt']]`. We use `map` for this.
* When the shell command is run `"cat {input} > {output}"`, the elements of the list are joined with a space. `cat input/A_1.txt input/A_2.txt > output/A_joined.txt`

### Other miscellaneous things
Directives like `params` and `threads` are fixed. Custom directives are not allowed.

Nested lists. `config["hey"]["oh"]`

Copying config file into output


## General guidelines
Link to the original Snakemake [tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html)

### Setting up
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

### Running on SLURM clusters (specific to me)
```
mkdir slurm/
snakemake -j 16 --cluster-config $CLST/cluster.json --cluster \
  "sbatch -J {cluster.jobname} -n {cluster.cores} -o {cluster.output} \
  -N {cluster.nodes}  -p {cluster.partition} \
  --mail-user {cluster.mailuser} --mail-type {cluster.mailtype}" 
```

### More advanced stuff
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