nextflow.enable.dsl=2

include { CELLTYPIST } from "$projectDir/workflows/celltypist"

samples_ch = Channel.fromPath(params.samples, checkIfExists: true)
    .splitCsv(header: true, sep: "\t", strip: true)
    .map( row -> [row.sample, row.countmtx_norm, row.cell_meta] )

workflow {
    println "Samplesheet: $params.samples"
    println "CellTypist model: $params.model"
    println "Output directory: $params.outdir\n"

    CELLTYPIST(samples_ch)
}
