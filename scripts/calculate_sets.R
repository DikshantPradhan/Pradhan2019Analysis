# GENERATE SETS

rm(list = ls())
# source('~/GitHub/PathwayMining/raptor_coupling.R')
# source('~/GitHub/PathwayMining/grb_tools.R')
# source('~/GitHub/PathwayMining/load_mod.R')
# source('~/GitHub/PathwayMining/gene_tools.R')
# source('~/GitHub/PathwayMining/set_tools.R')
# library(parallel)
# library(tictoc)
library(devtools)
library(roxygen2)
library(pathwaymining)
library(sybil)
library(parallel)
library(tictoc)
library(gurobi)
library(dplyr)
library(pryr)
library(data.tree)

model_data_generation <- function(sybil_model, grb_model, model_name,
                                  r0 = TRUE, gr0 = TRUE, g0 = TRUE, r1 = TRUE, gr1 = TRUE, g1 = TRUE, gpr_save = FALSE, cores = 1,
                                  directory = ''){
  
  # setup
  falcon_model <- GRB_generate_falcon_model(sybil_model)
  vars <- grb_model$varnames
  n <- length(vars)
  
  falcon_vars <- falcon_model$varnames
  falcon_n <- length(falcon_vars)
  non_vars <- grep('conversion', falcon_vars)
  
  non_gene_assc_rxns <- which(sybil_model@genes == "")
  gene_indexes <- grep('Ex_a', falcon_vars)
  falcon_rxn_idxs <- c(non_gene_assc_rxns, gene_indexes)
  
  falcon_full_rxn_idxs <- 1:(non_vars[1]-1)
  
  falcon_gene_rxn_idxs <- grep('Ex_a_', falcon_vars)
  if (length(falcon_rxn_idxs) == 0){
    falcon_rxn_idxs <- 1:(falcon_gene_rxn_idxs[length(falcon_gene_rxn_idxs)])
  }
  
  gpr <- generate_gpr(sybil_model)
  if (gpr_save){
    save(gpr, file = paste(directory, model_name, '_gpr.RData', sep = ''))
  }
  
  # R0 sets
  if (r0){
    print('r0 sets...')
    print(model_name)
    r0_coupling_mtx <- flux_coupling_raptor(grb_model)$coupled
    r0_sets <- get_list_of_sets_from_mtx(r0_coupling_mtx) #get_list_of_sets(return_couples(r0_coupling_mtx))
    save(r0_sets, file = paste(directory, model_name, '_r0_sets.RData', sep = ''))
    
    # GR0 sets
    if (gr0){
      print('gr0 sets...')
      gr0_sets <- gene_set_from_rxn_set(gpr, r0_sets)
      save(gr0_sets, file = paste(directory, model_name, '_gr0_sets.RData', sep = ''))
    }
  }
  
  # G0 sets
  if (g0){
    print('g0 sets...')
    print(model_name)
    g0_coupling_mtx <- flux_coupling_raptor(falcon_model, reaction_indexes = falcon_rxn_idxs)$coupled
    save(g0_coupling_mtx, file = paste(directory, model_name, '_g0_coupling_mtx.RData', sep = ''))
    g0_sets <- get_list_of_sets_from_mtx(g0_coupling_mtx)
    save(g0_sets, file = paste(directory, model_name, '_g0_sets.RData', sep = ''))
  }
  
  # R1 Sets
  if (r1){
    print('r1 sets...')
    print(model_name)
    coupling_vector_list <- GRB_generate_set_lists_cluster(grb_model, 1:n, compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = cores)
    save(coupling_vector_list, file = paste(directory, model_name, '_r1_coupling_vector_list.RData', sep = ''))
    r1_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector_list, n, vars) #coupling_matrix_from_array(mutans_coupling_array)
    r1_sets <- get_list_of_sets_from_mtx(r1_matrix)
    save(r1_sets, file = paste(directory, model_name, '_r1_sets.RData', sep = ''))
    
    # GR1 Sets
    if (gr1){
      print('gr1 sets...')
      gr1_sets <- gene_set_from_rxn_set(gpr, r1_sets)
      
      save(gr1_sets, file = paste(directory, model_name, '_gr1_sets.RData', sep = ''))
    }
  }
  
  # G1 Sets
  if (g1){
    non_gene_assc_rxns <- which(sybil_model@genes == "")
    gene_indexes <- grep('Ex_a', falcon_vars)
    suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)
    
    print('g1 sets...')
    print(model_name)
    coupling_vector_list <- GRB_generate_set_lists_cluster(falcon_model, suppression_idxs = suppr_indexes, reaction_indexes = suppr_indexes,
                                                           compare_known_init_sets = TRUE, optimize_suppr=TRUE, cores = cores)
    save(coupling_vector_list, file = paste(directory, model_name, '_g1_coupling_vector_list.RData', sep = ''))
    g1_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector_list, falcon_n, falcon_vars)
    g1_sets <- get_list_of_sets_from_mtx(g1_matrix)
    save(g1_sets, file = paste(directory, model_name, '_g1_sets.RData', sep = ''))
  }
}

load('pao_model.RData')
pao <- GRB_pao_model()
load('yeast_model.RData')
yeast <- GRB_yeast_model()
load('mutans_model.RData')
mutans <- GRB_mutans_model()
model_data_generation(pao_model, pao, 'pao', r0 = TRUE, gr0 = TRUE, g0 = TRUE, r1 = FALSE, gr1 = FALSE, g1 = FALSE)
model_data_generation(mutans_model, mutans, 'mutans', r0 = TRUE, gr0 = FALSE, g0 = TRUE, r1 = FALSE, gr1 = FALSE, g1 = FALSE)
model_data_generation(yeast_model, yeast, 'yeast', r0 = TRUE, gr0 = TRUE, g0 = TRUE, r1 = FALSE, gr1 = FALSE, g1 = FALSE)

# R1 Sets

# PAO
pao <- GRB_pao_model()
vars <- pao$varnames
n <- length(vars)

output_file <- 'pao_r1_coupling.csv'

coupling_vector <- read_coupling_csv(output_file)
avoid_rxns <- coupling_vector$completed_idxs

coupling_vector_list <- GRB_generate_set_lists_cluster(pao, compare_known_init_sets = TRUE, optimize_suppr=TRUE, 
                                                       optimize_rxns = TRUE, cores = 8, avoid_idxs = avoid_rxns, file_output = output_file)
save(coupling_vector_list, file = 'pao_r1_coupling_vector.RData')

# sets
output_file <- 'pao_r1_coupling.csv'
pao <- GRB_pao_model()
vars <- pao$varnames
n <- length(vars)
pao_r0_coupling_mtx <- flux_coupling_raptor(pao)$coupled
coupling_vector <- read_coupling_csv(output_file)
pao_r1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
pao_r1_sets <- get_list_of_sets_from_mtx(pao_r1_coupling_matrix, pao_r0_coupling_mtx)
save(pao_r1_sets, file = 'pao_r1_sets.RData')

# MUTANS
mutans <- GRB_mutans_model()
vars <- mutans$varnames
n <- length(vars)

output_file <- 'mutans_r1_coupling.csv'

coupling_vector <- read_coupling_csv(output_file)
avoid_rxns <- coupling_vector$completed_idxs

coupling_vector_list <- GRB_generate_set_lists_cluster(mutans, compare_known_init_sets = TRUE, optimize_suppr=TRUE, 
                                                       optimize_rxns = TRUE, cores = 8, avoid_idxs = avoid_rxns, file_output = output_file)
save(coupling_vector_list, file = 'mutans_r1_coupling_vector.RData')

# sets
output_file <- 'mutans_r1_coupling.csv'
mutans <- GRB_mutans_model()
vars <- mutans$varnames
n <- length(vars)
mutans_r0_coupling_mtx <- flux_coupling_raptor(mutans)$coupled
coupling_vector <- read_coupling_csv(output_file)
mutans_r1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
mutans_r1_sets <- get_list_of_sets_from_mtx(mutans_r1_coupling_matrix, mutans_r0_coupling_mtx)
save(mutans_r1_sets, file = 'mutans_r1_sets.RData')

# YEAST
yeast <- GRB_yeast_model()
vars <- yeast$varnames
n <- length(vars)

output_file <- 'yeast_r1_coupling.csv'

coupling_vector <- read_coupling_csv(output_file)
avoid_rxns <- coupling_vector$completed_idxs

coupling_vector_list <- GRB_generate_set_lists_cluster(yeast, compare_known_init_sets = TRUE, optimize_suppr=TRUE, 
                                                       optimize_rxns = TRUE, cores = 8, avoid_idxs = avoid_rxns, file_output = output_file)
save(coupling_vector_list, file = 'yeast_r1_coupling_vector.RData')

# sets
output_file <- 'yeast_r1_coupling.csv'
yeast <- GRB_yeast_model()
vars <- yeast$varnames
n <- length(vars)
yeast_r0_coupling_mtx <- flux_coupling_raptor(yeast)$coupled
coupling_vector <- read_coupling_csv(output_file)
yeast_r1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
yeast_r1_sets <- get_list_of_sets_from_mtx(yeast_r1_coupling_matrix, yeast_r0_coupling_mtx)
save(yeast_r1_sets, file = 'yeast_r1_sets.RData')

# G1 SETS

# YEAST MODEL
load('maranas_model_lipid_exch.RData')
yeast_falcon_model <- GRB_generate_falcon_model(yeast_model)

vars <- yeast_falcon_model$varnames
n <- length(vars)

non_gene_assc_rxns <- which(yeast_model@genes == "")
gene_indexes <- grep('Ex_a', yeast_falcon_model$varnames)
suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)

coupling_vector <- read_coupling_csv('yeast_g1_coupling.csv')
avoid_rxns <- coupling_vector$completed_idxs

coupling_vector_list <- GRB_generate_set_lists_cluster(yeast_falcon_model, suppression_idxs = suppr_indexes, 
                                                       reaction_indexes = suppr_indexes, compare_known_init_sets = TRUE, 
                                                       optimize_suppr=TRUE, cores = 6, avoid_idxs = avoid_rxns, 
                                                       file_output = 'yeast_g1_coupling.csv')
save(coupling_vector_list, file = 'yeast_g1_coupling_vector.RData')

yeast_g1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
yeast_g1_sets <- get_list_of_sets_from_mtx(yeast_g1_coupling_matrix)
save(yeast_g1_sets, file = 'yeast_g1_sets.RData')


# PAO MODEL

load('~/GitHub/PathwayMining/data/pao_model/pao_model.RData')

output_file <- 'pao_g1_coupling.csv'

pao_falcon_model <- GRB_generate_falcon_model(pao_model)

vars <- pao_falcon_model$varnames
n <- length(vars)

non_gene_assc_rxns <- which(pao_model@genes == "")
gene_indexes <- grep('Ex_a', pao_falcon_model$varnames)
suppr_indexes <- c(non_gene_assc_rxns, gene_indexes)

#coupling_vector <- read_coupling_csv(output_file)
avoid_rxns <- c() #coupling_vector$completed_idxs

coupling_vector_list <- GRB_generate_set_lists_cluster(pao_falcon_model, suppression_idxs = suppr_indexes, 
                                                       reaction_indexes = suppr_indexes, compare_known_init_sets = TRUE, 
                                                       optimize_suppr=TRUE, cores = 6, avoid_idxs = avoid_rxns, file_output = output_file)
save(coupling_vector_list, file = 'pao_coupling_vector.RData')

pao_g1_coupling_matrix <- coupling_matrix_from_coupling_vector_list(coupling_vector$coupling_vector, n, vars)
pao_g1_sets <- get_list_of_sets_from_mtx(pao_g1_coupling_matrix)
save(pao_g1_sets, file = 'pao_g1_sets.RData')
