#!/usr/bin/env nextflow

// ************************************************************************
// ********** MOFIFY HERE TO CHANGE INPUT AND OUTPUT DIRECTORIES **********
// ************************************************************************

params.inputDir = '/PATH/TO/input'
inputDirMAGs = params.inputDir + "/MAGs"
inputDirMetagenomes = params.inputDir + "/metagenomes"
params.outputDir = '/PATH/TO/output'
params.krakendbDir = '/PATH/TO/tmpfs'

params.threads = Runtime.runtime.availableProcessors() // maximum number available
params.memoryMapping = true
params.pathToNodesDmp = "/PATH/TO/nodes.dmp"

/* Available tools for gene prediction:
   "Augustus", "FragGeneScan", "GeneMarkS", "GeneMarkS2", "GeneMarkST", "GlimmerHMM", "MetaGeneAnnotator", 
   "MetaGeneMark", "MetaGeneMark2", "Phanotate", "Prodigal", "Pyrodigal", "Pyrodigalstat", "Snap"
*/
params.tools_for_org_group = """{
    "Archaea":["FragGeneScan", "GeneMarkST","MetaGeneMark"], 
    "Bacteria":["MetaGeneAnnotator","MetaGeneMark","Pyrodigal"], 
    "Eukaryota":["Augustus","Pyrodigal","Snap"], 
    "Viruses":["MetaGeneAnnotator","GeneMarkS","GeneMarkS2"]
    }"""


// Leaving these paths blank will assume that they have been added to PATH
// The lisence of GeneMark* software does not allow them to be packaged with InvestiGut
pathToGeneMarkS = ""     // "/PATH/TO/gmsn.pl"
pathToGeneMarkS2 = ""    // "/PATH/TO/gms2.pl"
pathToGeneMarkST = ""    // "/PATH/TO/gmst.pl"
pathToMetaGeneMark = ""  // "/PATH/TO/gmhmmp"
pathToMetaGeneMark2 = "" // "/PATH/TO/run_mgm.pl"

// ************************************************************************
// ************************* END OF USER SETTINGS *************************
// ************************************************************************

kraken2PredictDir = params.outputDir + "/1_kraken2predict"
sortedMagsMetagenomesDir = params.outputDir + "/2_sortedmagmetagenomes"
predictGenesDir = params.outputDir + "/3_predictgenes"
gffToFastaDir = params.outputDir + "/4_gfftofasta"
diamondDbDir = params.outputDir + "/5_diamonddb"

// PATH/TO/aa/bb/cc/INPUT_DIR/xx/yy/zz/FILE.fa.gz -> /xx/yy/zz/
params.getPathToFile = { path ->
    def relativePath = path.toString().substring(params.inputDir.length())
    def lastIndex = relativePath.lastIndexOf('/')
    return relativePath.take(lastIndex + 1)
}

// PATH/TO/FILE.fa.gz -> FILE
params.getFileNoExtension = { path ->
    def fullPath = path.toString()
    def filename = new File(fullPath).name
    return filename.replaceFirst(/\.(fa|fna|fasta)(\.gz)?$/, '')
}

process KRAKEN2_PREDICT {
    
    input:
    val metagenomesChannel
    
    publishDir "${kraken2PredictDir}/${params.getPathToFile(metagenomesChannel)}", mode: 'symlink'
    
    output:
    path "${metagenomesChannel.baseName}.krakentree", emit: krakentree_file
    path "${metagenomesChannel.baseName}.krakentaxid", emit: krakentaxid_file
    val metagenomesChannel, emit: metagenomesChannel

    script:
    """
    kraken2 \
        --db ${params.krakendbDir} \
        ${params.memoryMapping ? '--memory-mapping' : ''} \
        --report ${metagenomesChannel.baseName}.krakentree \
        --output ${metagenomesChannel.baseName}.krakentaxid \
        --use-names \
        --confidence 0.15 \
        --threads ${params.threads} \
        $metagenomesChannel
    """
    
}

process GET_CODON_TABLES {
    publishDir "${baseDir}/scripts", mode: 'symlink'

    input:
    path krakentaxid

    script:
    """
    python3 ${baseDir}/scripts/get_codon_tables.py ${baseDir} ${params.pathToNodesDmp} ${krakentaxid}
    """

    output:
    path "tax_table_pairs.pkl"
    val "", emit: throwaway

}

process SORT_METAGENOMES {

    input:
    path krakentree
    path krakentaxid
    val metagenomesChannel
    val done
    
    publishDir "${sortedMagsMetagenomesDir}/${params.getPathToFile(metagenomesChannel)}", mode: 'symlink'
    
    output: 
    tuple val(metagenomesChannel), path("*${params.getFileNoExtension(metagenomesChannel)}.{fa,fna,fasta}{,.gz}")
    
    script:
    """
    python3 ${baseDir}/scripts/sort_metagenomes.py $metagenomesChannel '${krakentree}' '${krakentaxid}' '${baseDir}'
    """
}

process SORT_MAGS {

    input:
    val magsChannel
    val done
    
    publishDir "${sortedMagsMetagenomesDir}/${params.getPathToFile(magsChannel)}", mode: 'symlink'
    
    output:
    tuple val(magsChannel), path("*${params.getFileNoExtension(magsChannel)}.{fa,fna,fasta}{,.gz}")


    script:
    """
    python3 ${baseDir}/scripts/sort_mags.py $magsChannel '${baseDir}'
    """
}

process PREDICT_GENES {

    input:
    tuple val(magsAndMetagenomesChannel), path(sorted_mags_and_metagenomes)
    
    publishDir "${predictGenesDir}/${params.getPathToFile(magsAndMetagenomesChannel)}", mode: 'symlink'
    
    output:
    tuple val(magsAndMetagenomesChannel), path("*.{gff3,gff}"), path(sorted_mags_and_metagenomes)

    script:
    """
    python3 ${baseDir}/scripts/predict_genes.py $magsAndMetagenomesChannel '${baseDir}' '${pathToGeneMarkS}' '${pathToGeneMarkS2}' '${pathToGeneMarkST}' '${pathToMetaGeneMark}' '${pathToMetaGeneMark2}' '${sorted_mags_and_metagenomes}' '${params.tools_for_org_group}'
    """
}

process GFF_TO_FASTA {

    input:
    tuple val(magsAndMetagenomesChannel), path(predictions), path(sorted_mags_and_metagenomes)
    
    publishDir "${gffToFastaDir}/${params.getPathToFile(magsAndMetagenomesChannel)}", mode: 'symlink'
    
    output:
    path "*.{fa,fna,fasta}"

    script:
    """
    python3 ${baseDir}/scripts/gff_to_fasta.py $magsAndMetagenomesChannel '${baseDir}' '${sorted_mags_and_metagenomes}' '${predictions}'
    """
}

process DIAMOND_DB {

    input:
    path fasta_files

    publishDir "${diamondDbDir}", mode: 'symlink'

    output:
    path "my_db.dmnd"

    script:
    """
    diamond makedb --in ${fasta_files} --db my_db
    """
}

workflow {
    metagenomesChannel = Channel.fromPath(inputDirMetagenomes+'/**.{fa,fna,fasta}{,.gz}')
    magsChannel = Channel.fromPath(inputDirMAGs+'/**.{fa,fna,fasta}{,.gz}')
    
    pred = KRAKEN2_PREDICT(metagenomesChannel).krakentaxid_file
    gettab = GET_CODON_TABLES(pred.collect()).throwaway
    results_ch = SORT_METAGENOMES(KRAKEN2_PREDICT.out.krakentree_file, KRAKEN2_PREDICT.out.krakentaxid_file, KRAKEN2_PREDICT.out.metagenomesChannel,gettab)
    results_ch2 = SORT_MAGS(magsChannel, gettab)
    results_ch3 = PREDICT_GENES(results_ch.mix(results_ch2))
    results_ch4 = GFF_TO_FASTA(results_ch3)
    diamond_db = DIAMOND_DB(results_ch4.collectFile())
}