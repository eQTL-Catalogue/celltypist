process CELLTYPE_CONCORDANCE {
    container "quay.io/peepk/celltypist:v1.0.0"
    publishDir "$params.outdir/celltype_concordance" , mode: "copy"

    input:
    path(countmtxs)
    path(cell_metas)

    script:
    """
    plot_cm.py "\$(echo $countmtxs)" "\$(echo $cell_metas)"
    """

    output:
    path("*.png")
    path("*.tsv")
}
