#!/usr/bin/env nextflow

// ************************************************************************
// ********** MOFIFY HERE TO CHANGE INPUT AND OUTPUT DIRECTORIES **********
// ************************************************************************

inputDir = '/home/horneflablinux/Documents/workflows/investigut_pipeline/input'
inputDirMAGs = inputDir + "/MAGs"
inputDirMetagenomes = inputDir + "/metagenomes"
outputDir = '/home/horneflablinux/Documents/workflows/investigut_pipeline/output'
krakendbDir = '/home/horneflablinux/Documents/workflows/investigut_pipeline/tmpfs'

params.threads = 14
params.memory_mapping = true

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

pathToGeneMarkS = "/home/horneflablinux/Documents/workflows/investigut_pipeline/pipeline/gene_prediction/genemark_suite_linux_64/gmsuite/gmsn.pl"
pathToGeneMarkS2 = "/home/horneflablinux/Documents/workflows/investigut_pipeline/pipeline/gene_prediction/gms2_linux_64/gms2.pl"
pathToGeneMarkST = "/home/horneflablinux/Documents/workflows/investigut_pipeline/pipeline/gene_prediction/gmst_linux_64/gmst.pl"
pathToMetaGeneMark = "/home/horneflablinux/Documents/workflows/investigut_pipeline/pipeline/gene_prediction/MetaGeneMark_linux_64/mgm/gmhmmp"
pathToMetaGeneMark2 = "/PATH/TO/MetaGeneMark2"

// ************************************************************************
// ************************* END OF USER SETTINGS *************************
// ************************************************************************

kraken2PredictDir = outputDir + "/1_kraken2predict"
sortedMagsMetagenomesDir = outputDir + "/2_sortedmagmetagenomes"
predictGenesDir = outputDir + "/3_predictgenes"
gffToFastaDir = outputDir + "/4_gfftofasta"
diamondDbDir = outputDir + "/5_diamonddb"

// PATH/TO/aa/bb/cc/INPUT_DIR/xx/yy/zz/FILE.fa.gz -> /xx/yy/zz/
params.getPathToFile = { path ->
    def relativePath = path.toString().substring(inputDir.length())
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
        --db ${krakendbDir} \
        ${params.memory_mapping ? '--memory-mapping' : ''} \
        --report ${metagenomesChannel.baseName}.krakentree \
        --output ${metagenomesChannel.baseName}.krakentaxid \
        --use-names \
        --confidence 0.15 \
        --threads ${params.threads} \
        $metagenomesChannel
    """
    
}

process GET_CODON_TABLES_OLD {
    publishDir "${baseDir}", mode: 'symlink'

    input:
    val throwaway

    script:
    """
    python3 ${baseDir}/get_codon_tables.py
    """

    output:
    path "tax_table_pairs.pkl"
    val "", emit: throwaway

}

process GET_CODON_TABLES {
    publishDir "${baseDir}", mode: 'symlink'

    input:
    path krakentaxid

    script:
    """
    python3 ${baseDir}/get_codon_tables.py ${baseDir} ${krakentaxid}
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
    python3 ${baseDir}/sort_metagenomes.py $metagenomesChannel '${krakentree}' '${krakentaxid}' '${baseDir}'
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
    python3 ${baseDir}/sort_mags.py $magsChannel '${baseDir}'
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
    python3 ${baseDir}/predict_genes.py $magsAndMetagenomesChannel '${baseDir}' '${pathToGeneMarkS}' '${pathToGeneMarkS2}' '${pathToGeneMarkST}' '${pathToMetaGeneMark}' '${pathToMetaGeneMark2}' '${sorted_mags_and_metagenomes}' '${params.tools_for_org_group}'
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
    python3 ${baseDir}/gff_to_fasta.py $magsAndMetagenomesChannel '${baseDir}' '${sorted_mags_and_metagenomes}' '${predictions}'
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