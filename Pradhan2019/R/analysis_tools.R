
#' create a data frame for making graphs 
#' @param model sybil model
#' @param gr0 g(r) sets
#' @param gr1 g(r*) sets
#' @return data frame
#' @seealso 
#' @export
#' @examples
#' 
gene_data_frame <- function(model, gr0, gr1){
  genes <- model@allGenes
  
  gpr_idxs <- lapply(genes, function(x){get_set_idx(x, model@genes)})
  gpr_promiscuity <- vapply(gpr_idxs, function(x){length(x)}, c(1))
  gr0_idxs <- lapply(genes, function(x){get_set_idx(x, gr0)})
  gr0_promiscuity <- vapply(gr0_idxs, function(x){length(x)}, c(1))
  gr1_idxs <- lapply(genes, function(x){get_set_idx(x, gr1)})
  gr1_promiscuity <- vapply(gr1_idxs, function(x){length(x)}, c(1))
  
  data.frame(genes = genes, gpr_promiscuity = gpr_promiscuity, gr0_promiscuity = gr0_promiscuity, gr1_promiscuity = gr1_promiscuity)
}

# CORRELATION ANALYSIS

#' clean names of genes in sets for easier analysis
#' @param sets sets
#' @param pattern character pattern to clean
#' @return sets with cleaned names
#' @seealso 
#' @export
#' @examples
#' 
clean_sets <- function(sets, pattern="PA\\d{4}") {
  genes_only <- function(x) {
    x %>%
      stringr::str_extract(pattern) %>%
      purrr::discard(is.na)
  }
  sets %>%
    purrr::map(genes_only) %>%
    purrr::discard(~length(.) == 0)
}

no_diag <- function(X) {
  diag(X) = NA
  x <- as.numeric(X)
  x[!is.na(x)]
}

#' plot correlation density of sets
#' @param C correlation matric
#' @param sets sets to analyze correlation of
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
plot_in_out <- function(sets, C = C) {
  plot_title <- dplyr::quo_name(dplyr::enquo(sets))
  all_genes <- dimnames(C)[[1]]
  set_genes <- intersect(unique(purrr::flatten_chr(sets)), all_genes)
  sets <- purrr::map(sets, intersect, set_genes)
  noset_genes <- setdiff(all_genes, set_genes)
  setC = C[set_genes, set_genes]
  
  nonsing_sets <- sets %>% purrr::discard(~length(.) < 2)
  in_corrs = nonsing_sets %>%
    purrr::map(~ no_diag(setC[.,.])) %>%
    purrr::flatten_dbl()
  out_corrs = nonsing_sets %>%
    purrr::map(~ as.numeric(setC[.,setdiff(set_genes, .)])) %>%
    purrr::flatten_dbl()
  
  plot(density(no_diag(C)), main=plot_title)
  lines(density(no_diag(setC)), col="red")
  lines(density(in_corrs), col="blue")
  lines(density(out_corrs), col="green")
}

# HEATMAPS

#' calculate average correlation values between pairs in each set
#' @param C correlation matrix
#' @param set sets to extract pairs from for correlations
#' @return list of integers
#' @seealso 
#' @export
#' @examples
#' 
set_cor <- function(set, avg = TRUE, C = C){
  if (length(set) < 2){
    return(0)
  }
  
  cor_list <- c()
  
  for (i in 1:(length(set)-1)){
    for (j in (i+1):length(set)){
      # print(paste(set[i],set[j]))
      if (set[i] %in% rownames(C) & set[j] %in% rownames(C)){
        cor_list <- c(cor_list, C[set[i],set[j]])
      }
    }
  }
  
  if (!avg){
    return(cor_list)
  }
  
  return(mean(cor_list))
}

g0_gr0_compar <- function(g0_set, gr0_set, C){
  new_gr0 <- gr0_set[!gr0_set %in% g0_set]
  genes <- c(g0_set, new_gr0)
  
  return(C[genes,genes])
}

g0_gr0_sets <- function(g0, gr0, g0_idx){
  # print(g0[[g0_idx]])
  
  first_rxn <- g0[[g0_idx]][1]
  gr0_idx <- get_set_idx(first_rxn, gr0)[1]
  # print(gr0[[gr0_idx]])
  
  return(c(g0_idx, gr0_idx))
}

#' plot correlation values across all sets
#' @param df dataframe of correlation values
#' @param title string indicating the title of the plot
#' @param x string indicating x-axis label
#' @param y string indicating y-axis label
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
plot_correlation_and_quantiles <- function(df, title = '', x = '', y = ''){
  df <- df[rev(order(df$avg_set_cor)),]
  df <- df[which(df$avg_set_cor != 0),]
  len <- 1:nrow(df)
  plot(0, type = 'n', xlim = c(0, nrow(df)), ylim = c(-0.5, 1), , main = title, xlab = x, ylab = y)
  for (i in len){
    rect(i-0.5, df$quantiles[[i]][2], i-0.5, df$quantiles[[i]][4], col = 'blue', border = 'blue')
  }
  lines(len, df$avg_set_cor, type = 'l', col = 'red')
}

#' calculate average correlation values in all sets as well as quantiles
#' @param sets list of sets
#' @param C gene expression correlation matrix
#' @return dataframe of various analyses of the sets
#' @seealso 
#' @export
#' @examples
#' 
set_analysis <- function(sets, C){
  avg_set_cor <- vapply(1:length(sets), function(x){return(set_cor(sets[[x]], C = C))}, FUN.VALUE = c(1))
  avg_set_cor <- avg_set_cor
  set_quantile <- lapply(1:length(sets), function(x){return(quantile(set_cor(sets[[x]], avg = FALSE, C = C)))})
  set_analysis <- data.frame(avg_set_cor)
  set_analysis$quantiles <- set_quantile
  
  return(set_analysis)
}

#' draw a heatmap of correlation values between genes in a specified G set and a G(R) set containing it
#' @param C gene expression correlation matrix
#' @param g0_idx integer indicating the index of the G set of interest
#' @param gr0_idx integer indicating the index of the G(R) set of interest
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
draw_heatmap_g0 <- function(g0_idx, gr0_idx, g0, gr0, C){
  idxs <- c(g0_idx, gr0_idx) #g0_gr0_sets(g0_idx)
  compar <- g0_gr0_compar(g0[[g0_idx]], gr0[[gr0_idx]], C)
  cols <- rep('black', nrow(compar))
  cols[row.names(compar) %in% g0[[idxs[1]]]] <- 'red'
  col_list <- c('brown3', 'cadetblue3','chartreuse3','chocolate3','burlywood3','darkgoldenrod4','aquamarine3','darkorchid3','deeppink','lightgoldenrod4',
                'midnightblue','mistyrose3','lightsteelblue2','lightsalmon','yellowgreen','tan1','springgreen3','snow3','yellow3','slateblue','sienna3',
                'violetred3','seagreen4','royalblue2','palevioletred','paleturquoise3','seagreen3','plum4')
  row_cols <- rep('white', nrow(compar))
  col_cols <- rep('white', nrow(compar))
  
  g0_sets <- vapply(rownames(compar), function(x){get_set_idx(x, g0)}, c(1))
  
  g0_colors <- sample(col_list, size = length(unique(g0_sets)))
  
  colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
  my_palette <- colorRampPalette(c("darkorange2", "white", "darkmagenta"))(n = 299)
  
  for (set in 1:length(unique(g0_sets))){
    col_cols[which(g0_sets == unique(g0_sets)[set])] <- g0_colors[set]
  }
  heatmap.2(compar,Rowv=FALSE,Colv=FALSE, RowSideColors = col_cols, ColSideColors = col_cols, col=my_palette, 
            symbreaks = TRUE, symm = TRUE, dendrogram = 'none', trace = 'none', density.info = 'none', key = TRUE, keysize = 2,
            key.xlab = "Correlation Coefficient")#, lhei = c(4,8), lwid = c(2,4))
  # text(c(1.65,0.3), "G0 Set", pos = 2, offset = 5.5, side = 1)
}

# SET ANALYSIS AND BOOTSTRAPPING

#' generate a data frame for use in analysis from gene sets
#' @param set G sets or G(R) sets
#' @return data frame
#' @seealso 
#' @export
#' @examples
#' 
gene_set_dataframe <- function(sets){
  num_genes <- lapply(sets, length)
  
  frame <- cbind(sets, sets, num_genes)
  colnames(frame) <- c('sets', 'clean sets', '# genes')
  return(frame)
}

#' return a matrix indicating the index of the set at which each gene can be found
#' @param genes list of n genes of interest
#' @param set list of sets
#' @return n by 1 matrix of integers
#' @seealso 
#' @export
#' @examples
#' 
map_genes_to_sets <- function(genes, set){
  map <- matrix(data = 0, nrow = length(genes), ncol = 1)
  for (i in 1:length(map)){
    pos <- grep(genes[i], set)
    if (length(pos) < 1){next}
    map[i] <- pos
  }
  return(map)
}

#' convert sets to groups of binary integers depending on whether or not each element fits certain criteria
#' @param sets original sets
#' @param assignment list of elements which fit the criteria (these will be set to 1 and the remaining will be set to 0)
#' @return list of sets in which elements are either 1 or 0
#' @seealso 
#' @export
#' @examples
#' 
binary_set_assignment <- function(sets, assignment){
  for (i in 1:length(sets)){
    if (is.null(sets[[i]])){next()}
    for (j in 1:length(sets[[i]])){
      if (sets[[i]][j] %in% assignment){
        sets[[i]][j] <- 1
      }
      else {
        sets[[i]][j] <- 0
      }
    }
  }
  
  return(sets)
}

#' calculate the fraction of each binary set whose elements are 1
#' @param binary_sets list of n sets whose elements are 1 or 0
#' @return n by 1 matrix of doubles
#' @seealso 
#' @export
#' @examples
#' 
binary_set_histogram <- function(binary_sets){
  assignment_fraction <- matrix(data = -1, nrow = length(binary_sets), ncol = 1)
  for (i in 1:length(binary_sets)){
    assignment_fraction[i] <- length(which(binary_sets[[i]] == 1))/length(binary_sets[[i]])
  }
  return(assignment_fraction)
}

double_histogram <- function(data1, data2, title){
  dat1 = data.frame(x = as.numeric(data1), group="observed")
  dat2 = data.frame(x = as.numeric(data2), group="bootstrapped")
  dat = rbind(dat1, dat2)
  ggplot(dat, aes(x, fill=group, colour=group)) +
    geom_histogram(aes(y=..density..), breaks=seq(0,1,0.05), alpha=0.6, 
                   position="identity", lwd=0.2) +
    ggtitle(title)
}

#' expand a dataframe of gene sets to include binary sets indicating genes which show differential fitness, differential expression, either, purely one or the other, or both  
#' @param set_df dataframe of gene sets
#' @param sets gene sets
#' @param df_genes list of genes which show differential fitness (df) effects
#' @param de_genes list of genes which show differential expression (de)
#' @param interest_genes list of genes which show differential fitness or expression
#' @param pure_df_genes list of genes which show df but not de 
#' @param pure_de_genes list of genes which show de but not df
#' @param dfde_genes list of genes which show df and de effects
#' @return dataframe of sets and binary sets for all listed criteria
#' @seealso 
#' @export
#' @examples
#' 
gene_set_expanded_dataframe <- function(set_df, sets, df_genes, de_genes, interest_genes, pure_df_genes, pure_de_genes, dfde_genes){
  df <- data.frame(set_df)
  df$df <- binary_set_assignment(sets, df_genes)
  df$de <- binary_set_assignment(sets, de_genes)
  df$interest_genes <- binary_set_assignment(sets, interest_genes)
  df$pure_df <- binary_set_assignment(sets, pure_df_genes)
  df$pure_de <- binary_set_assignment(sets, pure_de_genes)
  df$dfde <- binary_set_assignment(sets, dfde_genes)
  df$df_frac <- binary_set_histogram(df$df)
  df$de_frac <- binary_set_histogram(df$de)
  df$pure_df_frac <- binary_set_histogram(df$pure_df)
  df$pure_de_frac <- binary_set_histogram(df$pure_de)
  df$dfde_frac <- binary_set_histogram(df$dfde)
  df$interest_frac <- binary_set_histogram(df$interest_genes)
  
  return(df)
}

#' plot boostrapped values to assess the spread of set completeness
#' @param mtx matrix of bootstrapped set-fraction values
#' @param xlimits range of x-axis for plot
#' @param ylimits range of y-axis for plot
#' @param start first column to plot (will plot from start to ncol). the first column is a list of zeros so the default value is 2
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
plot_density <- function(mtx, xlimits=c(0, 1), ylimits=c(0, 200), start = 2){
  n <- ncol(mtx)
  plot(1, type="n", xlab="", ylab="", xlim=xlimits, ylim=ylimits)
  for (i in start:n){
    freqs <- table(as.numeric(mtx[,i]))
    freqs_x <- as.numeric(names(freqs))
    lines(freqs_x, freqs, col = 'grey')
  }
}

#' plot boostrapped values to assess the spread of set completeness against observed values
#' @param mtx matrix of bootstrapped set-fraction values
#' @param xlimits range of x-axis for plot
#' @param ylimits range of y-axis for plot
#' @param start first column to plot (will plot from start to ncol). the first column is a list of zeros so the default value is 2
#' @param freqs table of observed values
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
plot_density_compar <- function(mtx, xlimits=c(0, 1), ylimits=c(0, 200), start = 2, freqs){
  plot_density(mtx, xlimits = xlimits, ylimits = ylimits, start = start)
  freqs_x <- as.numeric(names(freqs))
  lines(freqs_x, freqs, col = 'blue')
}

#' plot the number of sets which are fully complete or empty
#' @param observed list of observed set fractions
#' @param simulated matrix of bootstrapped sete fractions
#' @param start first column of 'simulated' to plot (will plot from start to ncol). the first column is a list of zeros so the default value is 2 
#' @param end last column of 'simulated' to plot. the first column is a list of zeros so the default value is 1001 to hvae 1000 bootstrapped values
#' @param title plot title
#' @param color color for histogram bars
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
plot_set_characteristic_hypothesis <- function(observed, simulated, start = 2, end = 1001, title = '', color = 'lightskyblue'){
  obs_frac <- length(which(observed == 0 | observed == 1))/length(observed)
  sim_fracs <- c()
  
  for (i in start:end){
    sim <- length(which(simulated[,i] == 0 | simulated[,i] == 1))/nrow(simulated)
    sim_fracs <- c(sim_fracs, sim)
  }
  sim_fracs <- data.frame(sim_fracs)
  # print(length(which(sim_fracs < obs_frac)))
  
  colnames(sim_fracs) <- 'frac'
  ggplot(sim_fracs, aes(x=frac)) + #ggtitle(title) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.01,
                   colour="black", fill=color) +
    geom_density() + scale_x_continuous(limits = c(0,1)) +
    xlab('Complete Set Fraction') + ylab('Density') +
    geom_vline(xintercept = obs_frac, colour = 'red', linetype="dashed") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
}


#' function for one boostrap iteration for the distribution of dfde genes in gene sets
#' @param sets list of genes sets
#' @param n_df number of df genes to randomly assign
#' @param n_de number of de genes to randomly assign
#' @return list of bootstrapped dfde genes
#' @seealso 
#' @export
#' @examples
#' 
sample_dfde_sets <- function(sets, n_df, n_de){
  
  n_sets <- length(sets)
  all_genes <- unique(unlist(sets))
  n_genes <- length(unique(unlist(sets)))
  
  sample_sets <- function(sets, n){
    # shuffle sets
    set_order <- sample(1:n_sets, size = n_sets)
    sets <- sets[set_order]
    
    assignment <- rep(FALSE, n_genes)
    names(assignment) <- all_genes
    
    # set assignments
    remaining <- n
    set_idx <- 1
    while (remaining > 0){
      genes <- sets[[set_idx]]
      
      if (length(genes) > remaining){
        genes <- sample(genes, size = remaining)
        assignment[genes] <- TRUE
        remaining_df <- 0
      }
      
      assignment[genes] <- TRUE
      remaining <- remaining - length(genes)
      set_idx <- set_idx+1
    }
    return(assignment)
  }
  
  df_assignment <- sample_sets(sets, n_df)
  de_assignment <- sample_sets(sets, n_de)
  dfde_genes <- names(df_assignment[which(df_assignment & de_assignment)])
  
  return(dfde_genes)
}

#' bootstrap and plot the distribution of dfde genes
#' @param sets original gene sets
#' @param num_df number of df genes
#' @param num_de number of de genes
#' @param num_dfde original number of dfde genes
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
plot_dfde_sampling <- function(sets, num_df, num_de, num_dfde){
  n = 1000
  
  num_sampled_dfde <- c()
  for (i in 1:n){
    dfde_genes <- sample_dfde_sets(sets, num_df, num_de)
    num_sampled_dfde <- c(num_sampled_dfde, length(dfde_genes))
  }
  
  qplot(num_sampled_dfde, geom = 'histogram', bins = 50, main = 'dFdE gene bootstrap') + 
    geom_vline(xintercept = num_dfde, colour = 'red') + xlab('number of dfde genes') + ylab('frequency')
}

#' shuffle df and de set fractions and observe overlap to assign dfde
#' @param df_frac_list list of set-fractions of df genes
#' @param de_frac_list list of set-frctions of de genes
#' @param set_size_list list of integers that indicate the number of elements in each set
#' @param fracs boolean value indicating whether or not to return the fraction or number of bootstrapped dfde genes in each set. If TRUE, will return fractions.
#' @return list of integers or doubles indicating the bootstrapped fraction or number of dfde genes in each set 
#' @seealso 
#' @export
#' @examples
#' 
bootstrap_dfde <- function(df_frac_list, de_frac_list, set_size_list, fracs = TRUE){
  if (length(df_frac_list) != length(de_frac_list)){return()}
  n_sets <- length(df_frac_list)
  
  dfde_sample <- matrix(data = 0, ncol = 1, nrow = n_sets)
  
  for (i in 1:n_sets){
    df_frac <- df_frac_list[i]
    de_frac <- de_frac_list[i]
    
    set_size <- set_size_list[[i]]
    
    df_sampled_idxs <- sample(1:set_size, size = floor(set_size*df_frac))
    de_sampled_idxs <- sample(1:set_size, size = floor(set_size*de_frac))
    if (fracs){
      dfde_sample[i] <- length(intersect(df_sampled_idxs, de_sampled_idxs))/set_size
    }
    else {
      dfde_sample[i] <- length(intersect(df_sampled_idxs, de_sampled_idxs))
    }
  }
  return(dfde_sample)
}

#' bootstrap the distribution of dfde genes in sets by only shuffling the assignments within sets and plot
#' @param df dataframe of gene set information
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
plot_dfde_bootstrap <- function(df){
  reduced_df <- df[which(df$X..genes > 1),]
  n_sets <- nrow(df)
  
  n = 1000
  dfde_vec <- matrix(data = 0, nrow = n_sets, ncol = 1)
  
  for (i in 1:n){
    df_frac_sampled <- sample(df$df_frac, size = n_sets)
    de_frac_sampled <- sample(df$de_frac, size = n_sets)
    dfde_sample <- bootstrap_dfde(df_frac_sampled, de_frac_sampled, df$X..genes, fracs = TRUE)
    dfde_vec <- cbind(dfde_vec, dfde_sample)
  }
  
  plot_density(dfde_vec, ylimits = c(0,110))
  freqs <- table(as.numeric(reduced_df$dfde_frac))
  freqs_x <- as.numeric(names(freqs))
  lines(freqs_x, freqs, col = 'red')
  title(main="dFdE set composition", xlab="% of genes with fitness defect & expression change", ylab="frequency")
  
}

#' count the number of edges between df, de, and dfde sets
#' @param coupling_mtx matrix of boolean values indicating edges between sets indicated by rows and columns
#' @param df_frac list of doubles indicating the fraction of df genes in each set
#' @param de_frac list of doubles indicating the fraction of de genes in each set
#' @param dfde_var integer indicating how dfde edges should be interpreted [1 = exactly one set is dfde, 2 = at least one set is dfde, 3 = neither may be dfde, 4 = neither set is dfde, 5 = both sets are dfde, anything else defaults to 1]
#' @return n by 3 matrix (n being the total number of edges). First column indicates edges between df sets, second column indicates edges between de sets, third column indicates edges between df and de
#' @seealso 
#' @export
#' @examples
#' 
g1_architecture_measurement_binary <- function(coupling_mtx, df_frac, de_frac, dfde_var){
  edge_list <- matrix(nrow = length(which(coupling_mtx)), ncol = 3)
  set_idxs <- as.numeric(rownames(coupling_mtx))
  
  edge_idx <- 1
  for (i in 1:nrow(coupling_mtx)){
    for (j in which(coupling_mtx[i,])){
      set_idx_1 <- set_idxs[i]
      set_idx_2 <- set_idxs[j]
      df_diff <- df_frac[set_idx_1]&df_frac[set_idx_2]
      de_diff <- de_frac[set_idx_1]&de_frac[set_idx_2]
      # at least one dfde
      euclidian_diff <- ((df_frac[set_idx_1]&de_frac[set_idx_1])&(df_frac[set_idx_2]&!de_frac[set_idx_2]))|((df_frac[set_idx_1]&de_frac[set_idx_1])&(!df_frac[set_idx_2]&de_frac[set_idx_2]))|((df_frac[set_idx_2]&de_frac[set_idx_2])&(df_frac[set_idx_1]&!de_frac[set_idx_1]))|((df_frac[set_idx_2]&de_frac[set_idx_2])&(!df_frac[set_idx_1]&de_frac[set_idx_1]))
      if (dfde_var == 2){ # at least one dfde
        euclidian_diff <- ((df_frac[set_idx_1]&de_frac[set_idx_1])&(df_frac[set_idx_2]|de_frac[set_idx_2]))|((df_frac[set_idx_2]&de_frac[set_idx_2])&(df_frac[set_idx_1]|de_frac[set_idx_1]))
      }
      if (dfde_var == 3){
        euclidian_diff <- (df_frac[set_idx_1]&de_frac[set_idx_2])|(df_frac[set_idx_2]&de_frac[set_idx_1])
      }
      if (dfde_var == 4){ # neither dfde
        euclidian_diff <- ((df_frac[set_idx_1]&!de_frac[set_idx_1])&(de_frac[set_idx_2]&!df_frac[set_idx_2]))|((df_frac[set_idx_2]&!de_frac[set_idx_2])&(de_frac[set_idx_1]&!df_frac[set_idx_1]))
      }
      if (dfde_var == 5){ # both dfde
        euclidian_diff <- (df_frac[set_idx_1]&de_frac[set_idx_2])&(df_frac[set_idx_2]&de_frac[set_idx_1])
      }
      #((df_frac[set_idx_1]&de_frac[set_idx_1])&(df_frac[set_idx_2]&!de_frac[set_idx_2]))|((df_frac[set_idx_1]&de_frac[set_idx_1])&(!df_frac[set_idx_2]&de_frac[set_idx_2]))|((df_frac[set_idx_2]&de_frac[set_idx_2])&(df_frac[set_idx_1]&!de_frac[set_idx_1]))|((df_frac[set_idx_2]&de_frac[set_idx_2])&(!df_frac[set_idx_1]&de_frac[set_idx_1]))
      #((df_frac[set_idx_1]&de_frac[set_idx_1])&(df_frac[set_idx_2]|de_frac[set_idx_2]))|((df_frac[set_idx_2]&de_frac[set_idx_2])&(df_frac[set_idx_1]|de_frac[set_idx_1]))
      #(df_frac[set_idx_1]&de_frac[set_idx_2])|(df_frac[set_idx_2]&de_frac[set_idx_1]) 
      #((df_frac[set_idx_1]&!de_frac[set_idx_1])&(de_frac[set_idx_2]&!df_frac[set_idx_2]))|((df_frac[set_idx_2]&!de_frac[set_idx_2])&(de_frac[set_idx_1]&!df_frac[set_idx_1])) #(df_frac[set_idx_1]&de_frac[set_idx_2])|(df_frac[set_idx_2]&de_frac[set_idx_1])
      #(df_frac[set_idx_1]&de_frac[set_idx_2])&(df_frac[set_idx_2]&de_frac[set_idx_1])
      
      edge_list[edge_idx,1] <- df_diff
      edge_list[edge_idx,2] <- de_diff
      edge_list[edge_idx,3] <- euclidian_diff
      
      edge_idx <- edge_idx + 1
    }
  }
  
  return(edge_list)
}

#' boostrap the distribution of sets in the established G* architecture
#' @param coupling_mtx matrix indicating coupling between G sets
#' @param df_frac vector indicating the fraction of df genes in each set
#' @param de_frac vector indicating the fraction of de genes in each set
#' @param n number of bootstrap iterations to sample. Default is 1000
#' @return matrix of boolean values in which each row indicates the edge which is present of absent and each column represents a bootstrapped iteration and sampling variant
#' @seealso 
#' @export
#' @examples
#' 
g1_architecture_bootstrap <- function(coupling_mtx, df_frac, de_frac, n = 1000, dfde_var = 5){
  sets <- rownames(coupling_mtx)
  
  boot_sets <- sample(sets, size = length(sets))
  rownames(coupling_mtx) <- boot_sets
  colnames(coupling_mtx) <- boot_sets
  
  edge_list <- g1_architecture_measurement_binary(coupling_mtx, df_frac, de_frac, dfde_var = dfde_var)
  
  for (i in 2:n){
    boot_sets <- sample(sets, size = length(sets))
    rownames(coupling_mtx) <- boot_sets
    colnames(coupling_mtx) <- boot_sets
    
    edge_list <- cbind(edge_list, g1_architecture_measurement_binary(coupling_mtx, df_frac, de_frac, dfde_var = dfde_var))
  }
  
  return(edge_list)
}

#' plot histogram of the distribution of edges
#' @param edge_lists list of vectors of boolean elements indicating which sets have edges between them
#' @param og_edges vector of boolean values indicating the observed edges
#' @param xlimits x-axis range
#' @param title title of plot
#' @param color color of bars
#' @return NA
#' @seealso 
#' @export
#' @examples
#' 
binary_edge_histogram <- function(edge_lists, og_edges, xlimits = c(0,1), title = '', color = 'lightskyblue'){
  
  bootstrapped <- c()
  for (i in 1:ncol(edge_lists)){
    bootstrapped <- c(bootstrapped, length(which(edge_lists[,i]))/nrow(edge_lists))
  }
  
  observed <- length(which(og_edges))/length(og_edges)
  bootstrapped <- data.frame(bootstrapped)
  # print(length(which(bootstrapped < observed)))
  
  colnames(bootstrapped) <- 'frac'
  ggplot(bootstrapped, aes(x=frac)) + ggtitle(title) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.01,
                   colour="black", fill=color) +
    geom_density() + scale_x_continuous(limits = c(0,1)) +
    geom_vline(xintercept = observed, colour = 'red', linetype="dashed") + xlab('Fraction') + ylab('Count') +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
}
