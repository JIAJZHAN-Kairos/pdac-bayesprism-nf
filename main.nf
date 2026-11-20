nextflow.enable.dsl=2

// One BayesPrism run per bulk RNA source (mixtures from different quantification
// pipelines are deconvolved separately against the same single-cell reference).
process BAYESPRISM {
    tag "${source}"
    publishDir { "${params.outdir}/${source}" }, mode: 'copy'
    input:
    tuple val(source), path(mixture)
    path reference
    output:
    path '*'
    script:
    """
    export R_LIBS_USER="\$PWD/Rlib"
    install_bayesprism.R "\$R_LIBS_USER" ${task.cpus} > install.log 2>&1
    run_bayesprism_source.R '${source}' '${mixture}' '${reference}' . ${task.cpus}
    """
}

workflow {
    if (!params.mixtures || !params.reference || !params.outdir)
        error 'Set --mixtures, --reference and --outdir'
    mixtures = channel.fromPath(params.mixtures, checkIfExists: true).map { f ->
        tuple(f.baseName.replaceFirst(/^mixture_/, ''), f)
    }
    BAYESPRISM(mixtures, file(params.reference, checkIfExists: true))
}
