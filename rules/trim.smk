def is_single_end(name, replic, lane):
    return pd.isnull(fastq_info.loc[(name, replic, lane), "fq2"])

def get_fastq(wildcards):
    return indir + fastq_info.loc[(wildcards.name, wildcards.replic, wildcards.lane), ["fq1", "fq2"]].dropna()

rule trimSE:
    input:
        get_fastq
    output:
        fq = trimdir + "{name}-rep{replic}-{lane}.fq.gz"
    log:
        logdir + "trimmomatic/{name}-rep{replic}-{lane}.log"
    params:
        settings = config["trimmomatic"]["settings"]
    threads: 8
    conda:
        "../envs/preprocessing.yml"
    shell:
        "trimmomatic SE -threads {threads} {input} {output} {params.settings} 2>{log}"

rule trimPE:
    input:
        get_fastq
    output:
        fq1 = trimdir + "{name}-rep{replic}-{lane}.1.fq.gz",
        fq2 = trimdir + "{name}-rep{replic}-{lane}.2.fq.gz",
        fq1_unpaired = trimdir + "{name}-rep{replic}-{lane}.1.unpaired.fq.gz",
        fq2_unpaired = trimdir + "{name}-rep{replic}-{lane}.2.unpaired.fq.gz"
    log:
        logdir + "trimmomatic/{name}-rep{replic}-{lane}.log"
    params:
        settings=config["trimmomatic"]["settings"]
    threads: 8
    conda:
        "../envs/preprocessing.yml"
    shell:
        "trimmomatic PE -threads {threads} {input} {output.fq1} {output.fq1_unpaired} {output.fq2} {output.fq2_unpaired} {params.settings} 2>{log}"



