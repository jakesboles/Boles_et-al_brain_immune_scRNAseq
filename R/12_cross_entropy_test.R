#loading required libraries and setting up workspace
#load libraries
{
  library(Seurat)
  library(ggplot2)
  library(sctransform)
  library(SeuratWrappers)
  library(scales)
  library(dplyr)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
  library(rstatix)
  library(stringr)
  library(ggpattern)
  library(Matrix)
  
  library( digest )
  require( dunn.test )
  library( flowCore )
  library( FlowSOM )
  library( ggridges )
  require( RANN )
  library( RColorBrewer )
  library( reshape2 )
  library( Rtsne )
  library( umap )
  
  library(NMF)
}
#setting up workspace and location variables
{
  loc_workspace <- getwd()
  
  dir.create(c(paste0(loc_workspace, "cross_entropy_test/00_data/"),
               paste0(loc_workspace, "cross_entropy_test/00_scripts/"),
               paste0(loc_workspace, "cross_entropy_test/01_seurat_results/"),
               paste0(loc_workspace, "cross_entropy_test/02_tsne_diff/")))
  
  loc_data <- paste0(loc_workspace, "cross_entropy_test/00_data/")
  loc_scripts <- paste0(loc_workspace, "cross_entropy_test/00_scripts/")
  
  loc_seurat_results <- paste0(loc_workspace, "cross_entropy_test/01_seurat_results/")
  loc_tsne_diff <- paste0(loc_workspace, "cross_entropy_test/02_tsne_diff/")
  
}
#defining functions
{
  # function to generate colours for given length
  gg_colour_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # function to set random seed depending on base number and string
  set.seed.here <- function( seed.base, seed.char ){
    seed.add <- strtoi( substr( digest( seed.char, "xxhash32" ), 2, 8 ), 16 )
    seed.new <- seed.base + seed.add    
    set.seed( seed.new )
    invisible( seed.new )
  }
  
  # source( paste0( loc_scripts, "ce_diff_test.r" ) )
  # calculates test on differences of cross-entropy
  ce.diff.test <- function( 
    cross.entropy, 
    event.partition, 
    partition.label, partition.color, partition.line.type, 
    base.test, base.dist, 
    dendrogram.order.weight, 
    result, cdf.figure, dendrogram.figure 
  )
  {
    partition <- levels( event.partition )
    partition.n <- length( partition )
    
    cross.entropy.split <<- split( cross.entropy, event.partition ) #REMOVE SECOND <
    
    if ( base.test == "ks" )
    {
      if ( partition.n == 2 )
      {
        ks.test.res <- ks.test( cross.entropy.split[[ 1 ]], 
                                cross.entropy.split[[ 2 ]] )
        
        test.res <- list( ks.single = ks.test.res )
      }
      else if ( partition.n > 2 )
      {
        comparison.n <- partition.n * ( partition.n - 1 ) / 2
        
        ks.pair <- vector( "list", comparison.n )
        comparison <- character( comparison.n )
        D.stat <- numeric( comparison.n )
        p.value <- numeric( comparison.n )
        
        k <- 1
        
        for ( i in 1 : ( partition.n - 1 ) )
          for ( j in (i+1) : partition.n )
          {
            ks.pair[[ k ]] <- ks.test( cross.entropy.split[[ i ]], 
                                       cross.entropy.split[[ j ]] )
            
            comparison[ k ] <- sprintf( "%s - %s", 
                                        partition.label[ partition[ i ] ], 
                                        partition.label[ partition[ j ] ] )
            
            D.stat[ k ] <- ks.pair[[ k ]]$statistic
            
            p.value[ k ] <- ks.pair[[ k ]]$p.value
            
            k <- k + 1
          }
        
        p.value.adj <- p.adjust( p.value, "holm" )
        
        test.res <- list( ks.multiple = list( ks.pair = ks.pair, 
                                              comparison = comparison, D.stat = D.stat, 
                                              p.value = p.value, p.value.adj = p.value.adj ) )
      }
      else
        stop( "no partitions for testing cross-entropy differences" )
    }
    else if ( base.test == "rank" )
    {
      if ( partition.n == 2 )
      {
        wilcox.test.res <- wilcox.test( cross.entropy.split[[ 1 ]], 
                                        cross.entropy.split[[ 2 ]] )
        
        test.res <- list( rank.single = wilcox.test.res )
      }
      else if ( partition.n > 2 )
      {
        kruskal.test.res <- kruskal.test( cross.entropy, 
                                          event.partition )
        
        dunn.test.res <- dunn.test( cross.entropy, event.partition, 
                                    method = "holm", alpha = fcs.ce.diff.test.alpha, altp = TRUE, 
                                    kw = FALSE, table = FALSE, list = TRUE )
        
        test.res <- list( rank.multiple = list( 
          kruskal = kruskal.test.res, dunn = dunn.test.res ) )
      }
      else
        stop( "no partitions for testing cross-entropy differences" )
    }
    else
      stop( "wrong base test for cross-entropy differences" )
    
    if ( ! is.null( result ) )
    {
      result.file <- file( result, "w" )
      sink( result.file )
      
      tr.name <- names( test.res )
      stopifnot( length( tr.name ) == 1 )
      
      if ( tr.name == "ks.single" )
      {
        cat( "\n**  Kolmogorov-Smirnov test\n")
        
        print( test.res$ks.single )
      }
      else if ( tr.name == "ks.multiple" )
      {
        cat( "\n**  Multiple Kolmogorov-Smirnov tests with Holm correction\n\n")
        
        comparison.width <- max( nchar( test.res$ks.multiple$comparison ) )
        
        for ( i in 1 : length( test.res$ks.multiple$comparison ) )
          cat( sprintf( "%-*s\t\tD = %g\t\tpv = %g\t\tadj-pv = %g\n", 
                        comparison.width, 
                        test.res$ks.multiple$comparison[ i ], 
                        test.res$ks.multiple$D.stat[ i ], 
                        test.res$ks.multiple$p.value[ i ], 
                        test.res$ks.multiple$p.value.adj[ i ] ) )
      }
      else if ( tr.name == "rank.single" )
      {
        cat( "\n**  Wilcoxon rank sum test\n")
        
        print( test.res$rank.single )
      }
      else if ( tr.name == "rank.multiple" )
      {
        cat( "\n**  Kruskal-Wallis rank sum test\n")
        
        print( test.res$rank.multiple$kruskal )
        
        cat( "\n**  Dunn post-hoc test with Holm correction\n\n")
        
        comparison.width <- max( nchar( 
          test.res$rank.multiple$dunn$comparison ) )
        
        for ( i in 1 : length( test.res$rank.multiple$dunn$comparisons ) )
          cat( sprintf( "%-*s\t\tZ = %g\t\tpv = %g\t\tadj-pv = %g\n", 
                        comparison.width, 
                        test.res$rank.multiple$dunn$comparisons[ i ], 
                        test.res$rank.multiple$dunn$Z[ i ], 
                        test.res$rank.multiple$dunn$altP[ i ], 
                        test.res$rank.multiple$dunn$altP.adjusted[ i ] ) )
      }
      else
      {
        sink()
        close( result.file )
        stop( "unknown test in ce-diff result" )
      }
      
      sink()
      close( result.file )
    }
    
    if ( ! is.null( cdf.figure ) )
    {
      if ( is.null( partition.label ) )
        partition.label = partition
      
      if ( is.null( partition.color ) )
        partition.color <- rainbow( partition.n )
      
      if ( is.null( partition.line.type ) )
        partition.line.type <- rep( 1, partition.n )
      
      png( filename = cdf.figure, width = fcs.ce.diff.figure.cdf.width, 
           height = fcs.ce.diff.figure.cdf.height )
      
      par( mar = c( 5.5, 6, 2, 1.5 ) )
      
      plot( ecdf( cross.entropy ), ylim = c( 0, 1 ), 
            xlab = "Cross-entropy", ylab = "CDF", main = "", 
            cex.lab = 3, cex.axis = 2.5, 
            col = fcs.ce.diff.figure.cdf.all.color, 
            lwd = fcs.ce.diff.figure.line.width - 1, do.points = FALSE )
      
      for ( pall in partition )
      {
        ces <- cross.entropy.split[[ pall ]]
        ces.n <- length( ces )
        
        if ( ces.n < fcs.ce.diff.figure.cdf.resolution ) {
          plot( ecdf( ces ), col = partition.color[ pall ], 
                lty = partition.line.type[ pall ], 
                lwd = fcs.ce.diff.figure.line.width, 
                do.points = FALSE, add = TRUE )
        }
        else {
          ecdf.x <- sort( ces )
          ecdf.y <- 1 : ces.n / ces.n
          
          lines( ecdf.x, ecdf.y, col = partition.color[ pall ], 
                 lty = partition.line.type[ pall ], 
                 lwd = fcs.ce.diff.figure.line.width )
        }
      }
      
      legend( "bottomright", 
              legend = c( fcs.ce.diff.figure.cdf.all.label, partition.label ), 
              col = c( fcs.ce.diff.figure.cdf.all.color, partition.color ), 
              lty = c(1, partition.line.type), 
              lwd = fcs.ce.diff.figure.line.width, 
              cex = fcs.ce.diff.figure.font.size )
      
      dev.off()
    }
    
    if ( ! is.null( dendrogram.figure ) && partition.n > 2 )
    {
      cross.entropy.dist <- matrix( 0, nrow = partition.n, 
                                    ncol = partition.n )
      
      for ( i in 1 : ( partition.n - 1 ) )
        for ( j in (i+1) : partition.n )
        {
          if ( base.dist == "ks" )
            cross.entropy.dist[ i, j ] <- ks.test( 
              cross.entropy.split[[ i ]], 
              cross.entropy.split[[ j ]] 
            )$statistic
          else if ( base.dist == "median" )
            cross.entropy.dist[ i, j ] <- abs( 
              median( cross.entropy.split[[ i ]] ) - 
                median( cross.entropy.split[[ j ]] )
            )
          else
            stop( "wrong base dist for cross-entropy differences" )
          
          cross.entropy.dist[ j, i ] <- cross.entropy.dist[ i, j ]
        }
      
      cross.entropy.hclust <- hclust( as.dist( cross.entropy.dist ) )
      
      if ( ! is.null( dendrogram.order.weight ) )
        cross.entropy.hclust <- as.hclust( reorder( 
          as.dendrogram( cross.entropy.hclust ), 
          dendrogram.order.weight, 
          agglo.FUN = mean
        ) )
      
      if ( is.null( partition.label ) )
        partition.label = partition
      
      png( filename = dendrogram.figure, 
           width = fcs.ce.diff.figure.dendrogram.width, 
           height = fcs.ce.diff.figure.dendrogram.height )
      
      par( mar = c( 5, 5.6, 4, 1.4 ) )
      
      plot( cross.entropy.hclust, 
            labels = partition.label, hang = -1, 
            xlab = "", ylab = "", main = "", sub = "", cex.axis = 3, 
            cex = fcs.ce.diff.figure.font.size )
      
      dev.off()
    }
    
    test.res
  }
  
  # source( paste0( loc_scripts, "ce_diff_test_tsne.r" ) )
  # calculates cross-entropy for tsne plots and calls ce.diff.test
  ce.diff.test.tsne <- function( 
    orig.data, tsne.data, 
    event.partition, 
    partition.label = NULL, partition.color = NULL, partition.line.type = NULL, 
    base.test = "ks", base.dist = "ks", 
    prob.sample.n = NULL, dendrogram.order.weight = NULL, 
    result = NULL, cdf.figure = NULL, dendrogram.figure = NULL 
  )
  {
    stopifnot( nrow( orig.data ) == nrow( tsne.data ) && 
                 nrow( orig.data ) == length( event.partition ) )
    
    data.n <- nrow( orig.data )
    
    if ( ! is.null( prob.sample.n ) && prob.sample.n < data.n )
      prob.sample.idx <<- sample( data.n, prob.sample.n )
    else
      prob.sample.idx <<- 1 : data.n
    
    if ( fcs.use.cached.results && 
         file.exists( fcs.ce.diff.tsne.cache.file.path ) )
    {
      cat( "Using cached results for probability\n" )
      
      load( fcs.ce.diff.tsne.cache.file.path )
    }
    else
    {
      cat( "Calculating probability\n" )
      
      # sampling here temporary, until optimizing dist( tsne.dat ) below
      orig.tsne.prob <<- calculate.probability.tsne( 
        orig.data[ prob.sample.idx, ], 
        tsne.data[ prob.sample.idx, ] 
      )
      
      save( orig.tsne.prob, file = fcs.ce.diff.tsne.cache.file.path )
    }
    
    cross.entropy.all <<- calculate.cross.entropy( orig.tsne.prob$orig, 
                                                   orig.tsne.prob$tsne )
    #REMOVE SECOND <
    
    event.partition.all <<- event.partition[ prob.sample.idx ]
    
    ce.diff.test( 
      cross.entropy.all, 
      event.partition.all, 
      partition.label, partition.color, partition.line.type, 
      base.test, base.dist, 
      dendrogram.order.weight, 
      result, cdf.figure, dendrogram.figure
    )
  }
  
  calculate.probability.tsne <- function( orig.dat, tsne.dat )
  {
    orig.dat.n <- nrow( orig.dat )
    tsne.dat.n <- nrow( tsne.dat )
    
    stopifnot( orig.dat.n == tsne.dat.n )
    
    # find nearest neighbors in original space and their distances
    
    orig.dat.nn2 <- nn2( normalize_input( orig.dat ), 
                         k = fcs.ce.diff.tsne.perplexity.factor * fcs.tsne.perplexity + 1 )
    
    orig.dat.self.idx <- sapply( 1 : orig.dat.n, function( ri ) {
      ri.idx <- which( orig.dat.nn2$nn.idx[ ri, ] == ri )
      ifelse( length( ri.idx ) == 1, ri.idx, NA )
    } )
    
    stopifnot( ! is.na( orig.dat.self.idx ) )
    
    orig.neigh <- t( sapply( 1 : orig.dat.n, function( ri ) 
      orig.dat.nn2$nn.idx[ ri, - orig.dat.self.idx[ ri ] ] ) )
    
    orig.dist2 <- t( sapply( 1 : orig.dat.n, function( ri ) 
      orig.dat.nn2$nn.dists[ ri, - orig.dat.self.idx[ ri ] ]^2 ) )
    
    # calculate probabilities associated to distances in original space
    
    orig.stdev <- apply( orig.dist2, 1, function( dd2 ) {
      tsne.perplexity.error <- function( ss, dd2 ) {
        p <- exp( - dd2 / (2*ss^2) )
        if ( sum( p ) < .Machine$double.eps )
          p <- 1
        p <- p / sum( p )
        p <- p[ p > 0 ]
        2^( - sum( p * log2( p ) ) ) - fcs.tsne.perplexity
      }
      
      dd2.min.idx <- 1
      dd2.ascen <- sort( dd2 )
      while( dd2.ascen[ dd2.min.idx ] == 0 )
        dd2.min.idx <- dd2.min.idx + 1
      ss.lower <- dd2.ascen[ dd2.min.idx ]
      
      dd2.max.idx <- 1
      dd2.descen <- sort( dd2, decreasing = TRUE )
      while( is.infinite( dd2.descen[ dd2.max.idx ] ) )
        dd2.max.idx <- dd2.max.idx + 1
      ss.upper <- dd2.descen[ dd2.max.idx ]
      
      while( tsne.perplexity.error( ss.upper, dd2 ) < 0 )
      {
        ss.lower <- ss.upper
        ss.upper <- 2 * ss.upper
      }
      
      while( tsne.perplexity.error( ss.lower, dd2 ) > 0 )
      {
        ss.upper <- ss.lower
        ss.lower <- ss.lower / 2
      }
      
      uniroot( tsne.perplexity.error, dd2, 
               interval = c( ss.lower, ss.upper ), 
               tol = ( ss.upper - ss.lower ) * .Machine$double.eps^0.25 )$root
    } )
    
    orig.prob <- t( sapply( 1 : orig.dat.n, function( i ) {
      p <- exp( - orig.dist2[ i, ] /  ( 2 * orig.stdev[ i ]^2 ) )
      p / sum( p )
    } ) )
    
    # symmetrize probabilities in original space
    
    for ( i in 1 : orig.dat.n )
      for ( j2 in 1 : length( orig.neigh[ i, ] ) )
      {
        j <- orig.neigh[ i, j2 ]
        
        i2 <- match( i, orig.neigh[ j, ] )
        
        if ( ! is.na( i2 ) )
        {
          if ( j > i )
          {
            sym.prob <- ( orig.prob[ i, j2 ] + orig.prob[ j, i2 ] ) / 2
            orig.prob[ i, j2 ] <- sym.prob
            orig.prob[ j, i2 ] <- sym.prob
          }
        }
        else
          orig.prob[ i, j2 ] <- orig.prob[ i, j2 ] / 2
      }
    
    orig.prob <- sweep( orig.prob, 1, rowSums( orig.prob ), "/" )
    
    # get distances in tsne space for closest neighbors in original space
    
    tsne.dist2 <- t( sapply( 1 : tsne.dat.n, function( i )
      sapply( orig.neigh[ i, ], function( j )
        sum( ( tsne.dat[ i, ] - tsne.dat[ j, ] )^2 )
      )
    ) )
    
    # calculate probabilities associated to distances in tsne representation
    
    tsne.prob.factor <- tsne.dat.n / 
      ( 2 * sum( 1 / ( 1 + dist( tsne.dat )^2 ) ) )
    
    tsne.prob <- t( apply( tsne.dist2, 1, function( dd2 ) 
      p <- tsne.prob.factor / ( 1 + dd2 )
    ) )
    
    list( orig = orig.prob, tsne = tsne.prob )
  }
  
  calculate.cross.entropy <- function( prim.prob, secd.prob )
  {
    prim.prob.n <- nrow( prim.prob )
    secd.prob.n <- nrow( secd.prob )
    
    prim.prob.m <- ncol( prim.prob )
    secd.prob.m <- ncol( secd.prob )
    
    stopifnot( prim.prob.n == secd.prob.n && prim.prob.m == secd.prob.m )
    
    sapply( 1 : prim.prob.n, function( i ) 
      - sum( prim.prob[ i, ] * log( secd.prob[ i, ] ) )
    )
  }
  
  # source( paste0( loc_scripts, "plot_all_dmrd_figures.r" ) )
  # plots all dimensionality reduction figures
  plot.all.dmrd.figures <- function( 
    redu.data, redu.data.max, 
    redu.figure.lims.factor, redu.figure.point.size, 
    redu.figure.dir, redu.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition 
  )
  {
    redu.lims <- redu.figure.lims.factor * redu.data.max * c( -1, 1 )
    
    the.dmrd.figure.width <- fcs.dmrd.figure.width + 
      fcs.dmrd.label.factor.width * 
      max( nchar( fcs.condition.label ), nchar( flow.sample.label ), 
           nchar( fcs.cluster.label ) )
    
    the.dmrd.figure.width.multi <- 
      fcs.dmrd.figure.ncol * fcs.dmrd.figure.width + 
      fcs.dmrd.label.factor.width * 
      max( nchar( fcs.condition.label ), nchar( flow.sample.label ), 
           nchar( fcs.cluster.label ) )
    
    the.dmrd.figure.height <- fcs.dmrd.figure.height
    
    the.dmrd.figure.height.multi <- 
      fcs.dmrd.figure.nrow * fcs.dmrd.figure.height
    
    # plot all events colored by cluster
    
    redu.plot <- plot.dimensionality.reduction( 
      redu.data[ , 1 ], redu.data[ , 2 ], 
      event.partition = dmrd.event.cluster, 
      partition.label = fcs.cluster.label, 
      partition.color = fcs.cluster.color, 
      dmrd.lims = redu.lims, 
      point.size = redu.figure.point.size, 
      show.guide = TRUE 
    )
    
    ggsave( 
      file.path( redu.figure.dir, 
                 sprintf( "%s_all_events__cluster.png", redu.figure.plot ) ), 
      redu.plot, 
      width = the.dmrd.figure.width, height = the.dmrd.figure.height 
    )
    
    # plot all events colored by condition
    
    redu.plot <- plot.dimensionality.reduction( 
      redu.data[ , 1 ], redu.data[ , 2 ], 
      event.partition = dmrd.event.condition, 
      partition.label = fcs.condition.label, 
      partition.color = adjustcolor( fcs.condition.color, 
                                     alpha.f = fcs.dmrd.color.alpha ), 
      dmrd.lims = redu.lims, 
      point.size = redu.figure.point.size, 
      show.guide = TRUE 
    )
    
    ggsave( 
      file.path( redu.figure.dir, 
                 sprintf( "%s_all_events__condition.png", redu.figure.plot ) ), 
      redu.plot, 
      width = the.dmrd.figure.width, height = the.dmrd.figure.height 
    )
    
    # plot all events colored by each marker level
    
    for ( fch in fcs.channel )
    {
      dmrd.event.level <- dmrd.data[ , fch ]
      
      redu.plot <- plot.dimensionality.reduction( 
        redu.data[ , 1 ], redu.data[ , 2 ], 
        event.level = dmrd.event.level, 
        dmrd.lims = redu.lims, 
        point.size = redu.figure.point.size, 
        show.guide = TRUE, guide.name = fcs.channel.label[ fch ] 
      )
      
      ggsave( 
        file.path( redu.figure.dir, 
                   sprintf( "%s_all_events__%s.png", redu.figure.plot, 
                            fcs.channel.label[ fch ] ) ), 
        redu.plot, 
        width = the.dmrd.figure.width, height = the.dmrd.figure.height 
      )
    }
    
    # plot all conditions colored by cluster
    
    redu.plot <- plot.dimensionality.reduction( 
      redu.data[ , 1 ], redu.data[ , 2 ], 
      event.group = dmrd.event.condition, 
      group.label = fcs.condition.label, 
      event.partition = dmrd.event.cluster, 
      partition.label = fcs.cluster.label, 
      partition.color = fcs.cluster.color, 
      dmrd.lims = redu.lims, 
      dmrd.nrow = fcs.dmrd.figure.nrow, dmrd.ncol = fcs.dmrd.figure.ncol, 
      point.size = redu.figure.point.size, 
      show.guide = TRUE 
    )
    
    ggsave( 
      file.path( redu.figure.dir, 
                 sprintf( "%s_all_conditions__cluster.png", redu.figure.plot ) ), 
      redu.plot, 
      width = the.dmrd.figure.width.multi, 
      height = the.dmrd.figure.height.multi 
    )
    
    # plot all conditions colored by each marker level
    
    for ( fch in fcs.channel )
    {
      dmrd.event.level <- dmrd.data[ , fch ]
      
      redu.plot <- plot.dimensionality.reduction( 
        redu.data[ , 1 ], redu.data[ , 2 ], 
        event.group = dmrd.event.condition, 
        group.label = fcs.condition.label, 
        event.level = dmrd.event.level, 
        dmrd.lims = redu.lims, 
        dmrd.nrow = fcs.dmrd.figure.nrow, dmrd.ncol = fcs.dmrd.figure.ncol, 
        point.size = redu.figure.point.size, 
        show.guide = TRUE, guide.name = fcs.channel.label[ fch ] 
      )
      
      ggsave( 
        file.path( redu.figure.dir, 
                   sprintf( "%s_all_conditions__%s.png", redu.figure.plot, 
                            fcs.channel.label[ fch ] ) ), 
        redu.plot, 
        width = the.dmrd.figure.width.multi, 
        height = the.dmrd.figure.height.multi 
      )
    }
    
    # plot each condition colored by cluster
    
    for ( cond in fcs.condition )
    {
      dmrd.cond.idx <- which( dmrd.event.condition == cond )
      
      redu.cond.data <- redu.data[ dmrd.cond.idx, ]
      dmrd.cond.event.cluster <- dmrd.event.cluster[ dmrd.cond.idx ]
      
      redu.plot <- plot.dimensionality.reduction( 
        redu.cond.data[ , 1 ], redu.cond.data[ , 2 ], 
        event.partition = dmrd.cond.event.cluster, 
        partition.label = fcs.cluster.label, 
        partition.color = fcs.cluster.color, 
        dmrd.lims = redu.lims, 
        point.size = redu.figure.point.size, 
        show.guide = TRUE 
      )
      
      ggsave( 
        file.path( redu.figure.dir, 
                   sprintf( "%s_%s__cluster.png", redu.figure.plot, 
                            fcs.condition.label[ cond ] ) ), 
        redu.plot, 
        width = the.dmrd.figure.width, height = the.dmrd.figure.height 
      )
    }
    
    # plot each condition colored by sample
    
    for ( cond in fcs.condition )
    {
      dmrd.cond.idx <- which( dmrd.event.condition == cond )
      
      redu.cond.data <- redu.data[ dmrd.cond.idx, ]
      dmrd.cond.event.sample <- dmrd.event.sample[ dmrd.cond.idx ]
      
      redu.plot <- plot.dimensionality.reduction( 
        redu.cond.data[ , 1 ], redu.cond.data[ , 2 ], 
        event.partition = dmrd.cond.event.sample, 
        partition.label = flow.sample.label, 
        partition.color = adjustcolor( flow.sample.color.single, 
                                       alpha.f = fcs.dmrd.color.alpha ), 
        dmrd.lims = redu.lims, 
        point.size = redu.figure.point.size, 
        show.guide = TRUE 
      )
      
      ggsave( 
        file.path( redu.figure.dir, 
                   sprintf( "%s_%s__sample.png", redu.figure.plot, 
                            fcs.condition.label[ cond ] ) ), 
        redu.plot, 
        width = the.dmrd.figure.width, height = the.dmrd.figure.height 
      )
    }
  }
  
  
  # source( paste0( loc_scripts, "plot_dimensionality_reduction.r" ) )
  # plots one dimensionality reduction figure
  plot.dimensionality.reduction <- function( 
    dmrd.x, dmrd.y, 
    event.group = NULL, group.label = NULL, 
    event.partition = NULL, partition.label = NULL, partition.color = NULL, 
    event.level = NULL, 
    dmrd.lims = NULL, dmrd.nrow = NULL, dmrd.ncol = NULL, 
    point.size = NULL, show.guide = FALSE, guide.name = NULL 
  )
  {
    if ( ! is.null( event.partition ) )
    {
      if ( is.null( event.group ) )
        ggdf <- data.frame( dmrd.x, dmrd.y, event.partition )
      else
        ggdf <- data.frame( dmrd.x, dmrd.y, event.partition, event.group )
      
      if ( is.null( partition.label ) )
        plot.label <- waiver()
      else
        plot.label <- partition.label
      
      if ( show.guide )
        plot.guide <- guide_legend( keyheight = 0.8, 
                                    override.aes = list( size = fcs.dmrd.legend.point.size ), 
                                    label.theme = element_text( size = fcs.dmrd.legend.label.size ), 
                                    title = guide.name, title.position = "top", 
                                    title.theme = element_text( size = fcs.dmrd.legend.title.size ) )
      else
        plot.guide <- FALSE
      
      dmrd.plot <- ggplot( ggdf, aes( x = dmrd.x, 
                                      y = dmrd.y, color = event.partition ) ) + 
        scale_color_manual( values = partition.color, labels = plot.label, 
                            guide = plot.guide )
    }
    else if ( ! is.null( event.level ) )
    {
      if ( is.null( event.group ) )
        ggdf <- data.frame( dmrd.x, dmrd.y, event.level )
      else
        ggdf <- data.frame( dmrd.x, dmrd.y, event.level, event.group )
      
      if ( show.guide )
        plot.guide <- guide_colorbar( barwidth = 0.8, barheight = 10, 
                                      title = guide.name, title.position = "top", 
                                      title.theme = element_text( size = fcs.dmrd.legend.title.size ) )
      else
        plot.guide <- FALSE
      
      dmrd.plot <- ggplot( ggdf, aes( x = dmrd.x, y = dmrd.y, 
                                      color = event.level ) ) + 
        scale_color_gradientn( colors = fcs.dmrd.density.palette, 
                               labels = NULL, guide = plot.guide )
    }
    else
    {
      if ( is.null( event.group ) )
        ggdf <- data.frame( dmrd.x, dmrd.y )
      else
        ggdf <- data.frame( dmrd.x, dmrd.y, event.group )
      
      dmrd.plot <- ggplot( ggdf, aes( x = dmrd.x, y = dmrd.y ) )
    }
    
    if ( is.null( dmrd.lims ) ) {
      dmrd.xy.max <- max( abs( c( dmrd.x, dmrd.y ) ) )
      dmrd.lims <- dmrd.xy.max * c( -1, 1 )
    }
    
    if ( is.null( point.size ) )
      point.size <- 0.5
    
    dmrd.plot <- dmrd.plot + 
      coord_fixed() + 
      lims( x = dmrd.lims, y = dmrd.lims ) + 
      geom_point( shape = 20, stroke = 0, size = point.size ) + 
      theme_bw() + 
      theme( axis.title = element_blank(), 
             axis.text  = element_blank(), 
             axis.ticks = element_blank(), 
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank() )
    
    if ( ! is.null( event.group ) )
    {
      dmrd.plot <- dmrd.plot + 
        facet_wrap( vars( event.group ), 
                    labeller = as_labeller( group.label ), 
                    nrow = dmrd.nrow, ncol = dmrd.ncol ) + 
        theme( strip.background = element_rect( fill = "white" ), 
               strip.text = element_text( size = fcs.dmrd.group.title.size ) )
    }
    
    dmrd.plot
  }
  
}

#Start here----
##cross entropy tests
## for-looped through cell subsets
classes <- c("blymph", "cd45neg", "dendritic_cells", "granulocytes", 
             "macrophages", "microglia", "monocytes", "natural_killer", "tcells")
messages <- c("B-lymphocytes", "CD45- cells", "Dendritic cells", "Granulocytes",
              "Macrophages", "Microglia", "Monocytes", "Natural killer cells", "T-cells")
# calculate cross-entropy test for tsne by disease
for (i in seq_along(classes)){
  message(paste0("Running ", messages[i]))
  
  dir.create(paste0(loc_tsne_diff, classes[i]))
  
  obj <- readRDS(paste0("data_objects/08_", classes[i], ".RDS"))
  
  #testing it for samples/subjects/mice as the condition group
  fcs.condition <- as.character(sort(unique(obj@meta.data$Dissociation)))
  fcs.condition.n <- length(fcs.condition)
  fcs.condition.label <- fcs.condition
  names(fcs.condition.label) <- fcs.condition
  
  #some parameters copied as is from Carlos' script
  {
    # graphics parameters
    fcs.color.pool <- c( 
      brewer.pal( 8, "Set1" )[ -6 ], 
      brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
      adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
                   red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
      adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
                   red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
      adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
                   red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
      adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
                   red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
      adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
                   red.f = 0.3, green.f = 0.3, blue.f = 0.3 ), 
      adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
                   red.f = 0.3, green.f = 0.3, blue.f = 0.3 ) )
    fcs.color.pool.n <- length( fcs.color.pool )
    
    fcs.line.type.pool <- 1:6
    fcs.line.type.pool.n <- length( fcs.line.type.pool )
    
    fcs.condition.color <- rep( 
      fcs.color.pool, 
      ceiling( fcs.condition.n / fcs.color.pool.n ) 
    )[ 1 : fcs.condition.n ]
    names( fcs.condition.color ) <- fcs.condition
    
    fcs.condition.line.type <- rep( 
      fcs.line.type.pool, 
      ceiling( fcs.condition.n / fcs.line.type.pool.n ) 
    )[ 1 : fcs.condition.n ]
    names( fcs.condition.line.type ) <- fcs.condition
    
    fcs.ce.diff.figure.cdf.width <- 1200
    fcs.ce.diff.figure.cdf.height <- 800
    
    fcs.ce.diff.figure.dendrogram.width <- 1200
    fcs.ce.diff.figure.dendrogram.height <- 800
    
    fcs.ce.diff.figure.font.size <- 2
    fcs.ce.diff.figure.line.width <- 3
    
    fcs.ce.diff.figure.cdf.resolution <- 500
    fcs.ce.diff.figure.cdf.all.color <- "black"
    fcs.ce.diff.figure.cdf.all.label <- "All"
  }
  
  #global switches/parameters being used in the ce.diff.test.tsne function
  {
    fcs.use.cached.results <- TRUE
    fcs.ce.diff.tsne.perplexity.factor <- 3
    fcs.tsne.perplexity <- 30
    
    fcs.ce.diff.tsne.cache.file.path <- 
      paste0( loc_tsne_diff, classes[i], "/umap_all_group_cache.dat")
  }
  
  set.seed.here( 42, "calculate cross-entropy test for umap" )
  
  #only select PCs used for calculating t-SNE
  dmrd.data <- obj@reductions$pca@cell.embeddings[,1:30]
  tsne.data <- obj@reductions$umap@cell.embeddings
  
  dmrd.event.condition <- gsub(pattern = "_.*$", 
                               replacement = "", 
                               x = rownames(dmrd.data))
  
  gc()
  ce.diff.test.tsne.res.disease <- ce.diff.test.tsne(
    dmrd.data, tsne.data, 
    factor(x = obj@meta.data$Dissociation, 
           levels = unique(obj@meta.data$Dissociation)), 
    partition.label = fcs.condition.label, 
    partition.color = fcs.condition.color, # needs fcs.condition.n
    partition.line.type = fcs.condition.line.type, 
    base.test = "ks", #fcs.ce.diff.base.test, 
    base.dist = "ks", #fcs.ce.diff.base.dist, 
    prob.sample.n = NULL, #fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = c(1:fcs.condition.n), #fcs.ce.diff.figure.dendrogram.weight.condition, 
    result = paste0( loc_tsne_diff, classes[i], 
                     sprintf( "%s_group.txt", 
                              "/umap_ce_diff_result" #fcs.ce.diff.tsne.result 
                     ) ), 
    cdf.figure = paste0( loc_tsne_diff, classes[i], 
                         sprintf( "%s_group.png", 
                                  "/umap_ce_diff_cdf" #fcs.ce.diff.tsne.figure.cdf 
                         ) ), 
    dendrogram.figure = paste0( loc_tsne_diff, classes[i], 
                                sprintf( "%s_group.png", 
                                         "/umap_ce_diff_dendrogram" #fcs.ce.diff.tsne.figure.dendrogram 
                                ) ) 
  )
  
  text_info <- read.delim(file = paste0(loc_tsne_diff, classes[i],
                                        "/umap_ce_diff_result_group.txt"), 
                          stringsAsFactors = F, skip = 2)[2,1]
  
  DimPlot(object = obj, 
          reduction = "umap", 
          group.by = "Dissociation" 
  ) + 
    xlab(label = "Dimension 1") + ylab(label = "Dimension 2") + 
    labs(title = "", caption = text_info) + coord_fixed() + 
    theme(text = element_text(size = 20), 
          axis.text = element_text(size = 16))
  
}