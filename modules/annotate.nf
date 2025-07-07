process ANNOTATE {
    container "quay.io/peepk/celltypist:v1.0.0"
    publishDir "$params.outdir/annotate/$sample_id", pattern:"$countmtx_out" , mode: "copy"
    publishDir "$params.outdir/annotate/$sample_id", pattern:"*.tsv.gz" , mode: "copy"
    publishDir "$params.outdir/annotate/$sample_id", pattern:"*.png" , mode: "copy"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(countmtx), path(cell_meta)
    val model_name

    script:
    countmtx_out = "countmtx_w_celltypes.${sample_id}.h5ad"
    cell_meta_out = "cell_meta.${sample_id}.tsv"
    """
    annotate.py $countmtx $model_name $countmtx_out $cell_meta $sample_id $cell_meta_out
    mv figures/* .
    """

    output:
    path(countmtx_out), emit: annotated_countmtx
    path(cell_meta_out), emit: cell_meta
    path("*.tsv.gz")
    path("*.png")
}
