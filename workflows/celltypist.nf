include { ANNOTATE } from "$projectDir/modules/annotate"
include { CELLTYPE_CONCORDANCE } from "$projectDir/modules/concordance"

workflow CELLTYPIST {
    take:
    samples_ch

    main:
    // Add cell type labels to the count matrix with CellTypist
    ANNOTATE(samples_ch, params.model)
    annotated_countmtx_ch = ANNOTATE.out.annotated_countmtx.collect()
    cell_meta_ch = ANNOTATE.out.cell_meta.collect()

    // Draw a confusion matrix of cell type concordance with metadata
    CELLTYPE_CONCORDANCE(annotated_countmtx_ch, cell_meta_ch)
}
