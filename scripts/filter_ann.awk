# Filter a SnpEff-annotated VCF by the ANN INFO field.
#
# ANN subfields are: Allele|Annotation|Impact|Gene_Name|...
#
# Usage: awk -v mode=high|missense -f filter_ann.awk annotated.vcf
BEGIN { FS = "\t" }

/^#/ { print; next }

{
    keep = 0
    n = split($8, info, ";")
    for (i = 1; i <= n; i++) {
        if (info[i] !~ /^ANN=/) continue
        ann = substr(info[i], 5)
        m = split(ann, entries, ",")
        for (j = 1; j <= m; j++) {
            if (split(entries[j], f, "|") < 3) continue
            if (mode == "high" && f[3] == "HIGH") keep = 1
            if (mode == "missense" && f[2] ~ /missense/) keep = 1
        }
    }
    if (keep) print
}
