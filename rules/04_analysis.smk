rule ref_dict_fai:
    input:
        f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        fai=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
        dict=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
    conda:
        "../envs/04_analysis.yaml"
    shell:
        """
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input} -O {output.dict}
        """