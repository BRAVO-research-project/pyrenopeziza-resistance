

rm( list = ls() )
   
require('data.table', quietly = TRUE)
require('dplyr', quietly = TRUE)
require('tidyr', quietly = TRUE)


TEST = FALSE

VERSION = '1.1.0'

SRC.METAL = 'METAL'
SRC.TASSEL = 'TASSEL'

MODE.GEM = 'GEM'
MODE.SNP = 'SNP'

stopQuietly <- function(...) {
    blankMsg = sprintf( '\r%s\r', paste(rep(' ', getOption('width')-1L), collapse=' '));
    stop( simpleError( blankMsg ) );
}

getFileNames <- function() {
    if (src.type == SRC.METAL) {
        msg = paste('Select', mode, 'output file from', src.type)
        src.files = tcltk::tk_choose.files(caption = msg, multi = FALSE)
        if ( length( src.files ) == 0 ) stopQuietly()
    }
    if (src.type == SRC.TASSEL) {
        msg = paste('Select the *two*', mode, 'files from', src.type, '(use Ctrl+Click)')
        src.files = tcltk::tk_choose.files(caption = msg, multi = TRUE)
        if ( length( src.files ) == 0 ) stopQuietly()
        if ( length( src.files ) != 2 ) {
            stop( 'You must select 2 files, i.e. the effects file and the stats file from TASSEL!', call. = FALSE )
        }
    }

    msg = 'Select the Arabidopsis annotation and AC direction files'
    meta.files = tcltk::tk_choose.files(caption = msg, multi = TRUE)
    if ( length( meta.files ) == 0 ) stopQuietly()
    if ( length( meta.files ) != 2 ) {
        stop( 'You must select 2 files - the Arabidopsis annotation and the AC directions', call. = FALSE )
    }

    l = list( 'src.files' = src.files, 'meta.files' = meta.files )
    return(l)
}

getMetalData <- function( metal_file ) {
    # Load and clean the METAL data

    print( 'Loading METAL data...' )
    metal.dt = fread( metal_file )

    setnames( metal.dt, old = 'P-value', new = 'pvalue' )
    
    # remove bad p-values
    metal.dt = metal.dt[ !is.nan( pvalue ), ]
    
    # want log(p)
    metal.dt$log10P = -log10( metal.dt$pvalue )

    # SNPs will have ',X' at the end so remove it, if present
    metal.dt = separate(
        metal.dt,
        col = 'MarkerName',
        into = c('Marker'),
        sep = ',',
        extra = 'drop'
    )

    metal.dt = metal.dt[, c('Marker', 'Allele1', 'Allele2', 'log10P')]

    return(metal.dt)
}

getTasselData <- function( tassel_files ) {
    ## Get the output data from TASSEL, merge and append '-log10(p)'
    
    print( 'Loading TASSEL data...' )
    
    # 'tassel_files' has filenames for the 'stats' and 'effects' files
    # but they could be in any order so treat generally
    dt.1 = fread( tassel_files[1] )
    dt.2 = fread( tassel_files[2] )

    traits = sort( unique( dt.1$Trait ) )
    selected_trait = select.list(choices = traits,
                                 title = 'Traits',
                                 graphics = FALSE)

    new.prefix = paste0( prefix, '_', selected_trait )
    
    # add the trait to the output prefix
    assign("prefix", new.prefix, envir = .GlobalEnv)
    
    dt.1 = dt.1[ Trait == selected_trait ]
    dt.2 = dt.2[ Trait == selected_trait ]
    
    tassel.dt = merge( dt.1, dt.2, by = c('Trait', 'Marker') )
    tassel.dt = tassel.dt[ !is.nan( p ), ]
    tassel.dt$log10P = -log10( tassel.dt$p )

    tassel.dt$Allele1 = tassel.dt$Allele
    tassel.dt$Allele2 = ''

    if ( 'Estimate' %in% colnames( tassel.dt) ) {
        setnames( tassel.dt, old = 'Estimate', new = 'Effect' )
    }

    tassel.dt = tassel.dt[, c('Trait', 'Marker', 'Allele1', 'Allele2', 'log10P', 'Obs', 'Effect')]

    return(tassel.dt)
}

loadMetaData <- function() {
    print( 'Loading the annotation and the direction files...')
    
    dt.1 = fread( filenames$meta.files[1] )
    dt.2 = fread( filenames$meta.files[2] )

    annos.cols = c('Unigene', 'Unigene.Marker', 'Chr', 'AGI')
    graph.cols = c('Unigene', 'Unigene.Marker', 'Chr', 'sort.SNP', 'sort.GEM')
    
    if ( all( annos.cols %in% colnames( dt.1 ) ) ) {
        annos.dt = dt.1
        graph.dt = dt.2
    }

    if ( all( graph.cols %in% colnames( dt.1 ) ) ) {
        annos.dt = dt.2
        graph.dt = dt.1
    }

    old.name = ifelse( mode == MODE.SNP, 'Unigene.Marker', 'Unigene' )
    setnames( annos.dt, old = old.name, new = 'Marker' )
    setnames( graph.dt, old = old.name, new = 'Marker' )

    old.name = ifelse( mode == MODE.SNP, 'sort.SNP', 'sort.GEM' )
    setnames( graph.dt, old = old.name, new = 'sort' )
    
    annos.dt = annos.dt[, c('Marker', 'Chr', 'AGI')]
    annos.dt = annos.dt %>% drop_na(Marker)
    annos.dt = annos.dt[ !duplicated(annos.dt), ]
    
    graph.dt = graph.dt[, c('Marker', 'Chr', 'sort')]
    graph.dt = graph.dt %>% drop_na(Marker, sort)
    graph.dt = graph.dt[ !duplicated(graph.dt), ]
    
    meta.dt = merge( graph.dt, annos.dt, by = c('Marker', 'Chr'))

    chromo.dt = getChromoInfo( meta.dt )
    meta.dt = merge( meta.dt, chromo.dt, by = 'Chr' )
    
    return(meta.dt)
}

getChromoInfo <- function( graph_data ) {
    print( 'Summarising the chromosome information...')
    
    min_vals = setNames( aggregate( sort ~ Chr, data = graph_data, FUN = min ), nm = c( 'Chr', 'Start' ) )
    max_vals = setNames( aggregate( sort ~ Chr, data = graph_data, FUN = max ), nm = c( 'Chr', 'End' ) )
    
    chromo.dt = data.table( merge( min_vals, max_vals, by = 'Chr' ) )
    chromo.dt = separate( data = chromo.dt,
                          col = 'Chr',
                          into = c( 'Genome', 'ChromoNum' ),
                          sep = c( 1 ),
                          remove = FALSE,
                          convert = TRUE )
    
    return(chromo.dt)
}

loadData <- function() {
    meta.dt = loadMetaData()
    
    if ( src.type == SRC.METAL ) {
        all.dt = getMetalData( metal_file = filenames$src.files )
    } else {
        all.dt = getTasselData( tassel_files = filenames$src.files )
    }

    ## merge and sort...
    all.dt = merge( all.dt, meta.dt, by = 'Marker', na.rm = TRUE )
    setorder( all.dt, sort )

    return(all.dt)
}

createGenomeManhattan <- function(filename, data) {
    print( 'Generating Manhattan plot for genome...')
    
    ylim = c(0, max( data$log10P ) )

    # choose A and/or C, if present
    genome_labels = c()
    for (genome_label in c('A', 'C'))
        if ( genome_label %in% data$Genome )
            genome_labels = c( genome_labels, genome_label )
    
    jpeg(
        filename = filename,
        quality = 1000,
        width = 3000,
        height = 500 * length( genome_labels ),
        units = 'px',
        bg = 'white',
        res = NA
    )

    # Sets up a split screen - required when data is for A and C genomes
    par(mfrow = c( length(genome_labels), 1 ) )

    for (genome_label in genome_labels) {
        genome.dt = data[ Genome == genome_label, c('Marker', 'sort', 'log10P', 'ChromoNum', 'End')]
        genome.dt = genome.dt[ !duplicated(genome.dt), ]
        
        xlim = c(1, max(genome.dt$End))

        plot (
            genome.dt$sort,
            genome.dt$log10P,
            type = 'p',
            pch = 20,
            xlim = xlim,
            ylim = ylim,
            main = genome_label,
            cex.main = 2,
            ylab = '-log10(p)',
            xlab = ' ',
            cex.lab = 1.5,
            col = ifelse( genome.dt$ChromoNum %% 2 == 0, 'red', 'black' ),
            xaxt = 'n'
        )
    }

    close.screen(all = TRUE)
    dev.off()

    print( paste( 'Saved Manhattan plot to', filename ) )

    return()
}

createChromosomeManhattan <- function(data) {
    ## plot chosen graph and allow point identification

    print( 'Generating Manhattan plot for chromosomes. Please choose one from the list.')
    
    ylim = c(0, max(data$log10P))
    
    # get and sort the available chromosomes - want A10 to be last, if it's present
    choice.dt = data[, c('Chr', 'Genome', 'ChromoNum') ]
    choice.dt = choice.dt[ !duplicated(choice.dt), ]
    setorder( choice.dt, Genome, ChromoNum )
    choices = choice.dt$Chr
    
    while (TRUE) {
        choice = select.list(choices = choices,
                             title = 'Chromosome',
                             graphics = TRUE)

        if ( choice == '' ) {
            close.screen(all = TRUE)
            dev.off()
            break
        }

        chrom.dt = data[ Chr == choice ]
        plot.dt = chrom.dt[, c('Marker', 'sort', 'log10P', 'ChromoNum', 'Start', 'End')]
        plot.dt = plot.dt[ !duplicated(plot.dt), ]
        
        labels = plot.dt$Marker
        x = plot.dt$sort
        y = plot.dt$log10P

        # 'Start' and 'End' are all the same
        xlim = c( plot.dt$Start[1], plot.dt$End[1] )

        plot(
            x = x,
            y = y,
            type = 'p',
            pch = 20,
            xlim = xlim,
            ylim = ylim,
            main = choice,
            ylab = '-log10(p)',
            xlab = '',
            col = 'black',
            xaxt = 'n'
        )

        print( 'Select points on plot - press `escape` when finished')

        idx = identify(
            x = x,
            y = y,
            labels = labels,
            pos = TRUE,
            n = 10,
            plot = TRUE,
            atpen = TRUE,
            tolerance = 0.1,
            col = 'red',
            cex = 0.6
        )

        if (select.list( choices = c('Yes', 'No'), title = 'Save as jpeg?', graphics = TRUE ) == 'Yes') {
            filename = paste0 (prefix, '_', choice, '.jpeg')
            dev.copy(jpeg, filename)
            dev.off()

            print( paste('Plot saved to', filename ))
        }

        results = getHitsAndEffects( dt = chrom.dt, selectedMarkers = plot.dt[ idx$ind, 'Marker'] )
        print( results )

        csv_file = paste0( prefix, '_', choice, '.csv')
        write.csv( x = results, file = csv_file, row.names = FALSE )

        print( paste( 'Information on selected points saved to', csv_file ) )
    }
    
    return()
}

getHitsAndEffects <- function( dt, selectedMarkers ) {
    columns = c( 'Marker', 'Chr', 'AGI', 'Allele1', 'Allele2', 'log10P' )
    if ( src.type == SRC.TASSEL ) {
        columns = c( 'Trait', columns, 'Obs', 'Effect' )
    }
    
    dt.result = merge( dt, selectedMarkers, by = 'Marker' )
    dt.result = dt.result[, ..columns]
    dt.result = dt.result[ order( dt.result$log10P, decreasing = TRUE ), ]

    return(dt.result)
}

####################################################################################################
####################################################################################################

print( paste0( '*** Welcome to Grapher v', VERSION, ' ***') )
print( paste0( 'Working in ', getwd() ) )

src.type = select.list( choices = c(SRC.METAL, SRC.TASSEL), 
                        title = 'Please select the data source...', 
                        graphics = FALSE)
mode = select.list( choices = c(MODE.GEM, MODE.SNP), 
                    title = 'Please select the data type...', 
                    graphics = FALSE)
prefix = readline(prompt = 'Please enter output prefix...   ')

filenames = getFileNames()
source.dt = loadData()

createGenomeManhattan( filename = paste0( prefix, '_Manhattan.jpeg' ),
                       data = source.dt )

createChromosomeManhattan( data = source.dt )

print( paste0( '*** Goodbye ***') )
