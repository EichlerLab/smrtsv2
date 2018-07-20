rule plot_genotypes:
    input:
        'samples/{sample}/genotypes.tab'
    output:
        'samples/{sample}/genotypes.pdf'
    shell:
        """Rscript plot_genotypes.R {input} {output} {wildcards.sample}"""
