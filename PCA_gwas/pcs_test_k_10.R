setwd("/dscrhome/yy222/gas-rgls-master/scripts")
# for terminal options
library(optparse)

############
### ARGV ###
############

# define options
option_list = list(
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE, 
              help = "create plot of existing data"),
  make_option("--simple", action = "store_true", default = FALSE, 
              help = "simple plot version (for a grant proposal)"),
  make_option(c("-l", "--lite"), action = "store_true", default = FALSE, 
              help = "run lite version (focuses on best and least redundant methods)"),
  make_option(c("-n", "--n_ind"), type = "integer", default = 100, 
              help = "number of individuals", metavar = "int"),
  make_option(c("-m", "--m_loci"), type = "integer", default = 100000, 
              help = "number of loci", metavar = "int"),
  make_option(c("-k", "--k_subpops"), type = "integer", default = 10, 
              help = "admixture intermediate subpopulations", metavar = "int"),
  make_option(c("-f", "--fst"), type = "double", default = 0.3, 
              help = "FST (fixation index)", metavar = "double"),
  make_option(c("--bias_coeff"), type = "double", default = 0.5, 
              help = "admixture bias coeff", metavar = "double"),
  make_option(c("-g", "--generations"), type = "integer", default = 1, 
              help = "number of generations, for realistic local kinship", metavar = "int"),
  make_option("--herit", type = "double", default = 0.8, 
              help = "heritability", metavar = "double"),
  make_option("--m_causal", type = "integer", default = 100, 
              help = "num causal loci", metavar = "int"),
  make_option(c("-t", "--threads"), type = "integer", default = 1, 
              help = "number of threads (affects GCTA only)", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
generations <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
lite <- opt$lite
threads <- opt$threads
simple <- opt$simple

  library(readr) # to read style
  


  
  # load style table (maps methods to colors and line styles)
  style <- read_tsv('style.txt', col_types = 'cci')

  
  # else run simulation
  
  library(readr)       # to write kinship matrix
  library(tibble)      # to store data
  library(genio)       # to write BED files for external software
  library(popkin)      # to estimate kinship in RGLS
  library(popkinsuppl) # for PCA's kinship estimator
  library(lfa)         # GWAS gcatest
  library(gcatest)     # GWAS gcatest
  #library(qvalue)      # multiple hypothesis tests
  
  # standard code for a complex trait and an admixed population
  source('sim_geno_trait_k3.R')
  
  # load new functions from external scripts
  source('kinship_to_evd.R')
  source('gas_lm_optim.R')
  source('gas_pca_optim.R')
  source('gas_rgls.R')
  source('gas_lmm_gemma.R')
  source('gas_lmm_emmax.R')
  source('gas_lmm_gcta.R')
  source('paths.R')
  source("gas_plots.R")
  ############
  ### SIMS ###
  ############
     sequence<-seq(100,295,10)
  M_rmsd<-matrix(0,5,300)
  M_auc<-matrix(0,5,300)
  style <- read_tsv('style.txt', col_types = 'cci')
  for (pcs in sequence){
    for (j in 1:5){
      sim_geno_trait_k3 <- function(
        # params of population structure
        n_ind = 100, # number of individuals (columns)
        m_loci = 100000, # number of loci (rows)
        k_subpops = 10, # number of intermediate subpopulations
        fst = 0.3, # desired final Fst
        bias_coeff = 0.5, # bias coeff of standard Fst estimator
        # the family code gets triggered if generations > 1
        generations = 1, # number of generations (1 is unrelated founders)
        # params of complex trait
        m_causal = 100,
        herit = 0.8,
        verbose = TRUE, # to print messages
        low_mem = FALSE # for bnpsd::draw_all_admix
      ) {
        # this function is for generating random genotypes from an admixed population with 3 source subpopulations
        # meant to abstractly resemble Hispanics and African-Americans
        
        # things that can't be changed:
        # - Fst of intermediate subpopulations is a ramp proportional to 1:k
        # - Admixture proportions are from 1D geography model
        
        ####################
        ### 1: GENOTYPES ###
        ####################
        
        # simulating genotypes from an admixed population
        if (verbose) message('sim_geno_k3')
        
        # define population structure
        labs<-ceiling((1:1000)*(10/1000))
        # in this case return value is a named list with three items:
        admix_proportions <- admix_prop_indep_subpops(labs=labs) # admixture proportions
        inbr_subpops <- rep(fst,10) # rescaled Fst vector for intermediate subpops
        
        # draw allele freqs and genotypes
        if (verbose)
          message('draw_all_admix')
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, low_mem = low_mem)
        X <- out$X # genotypes
        p_anc <- out$p_anc # ancestral AFs
        
        # for the first pass of RGLS, let's calculate true kinship matrix
        if (verbose)
          message('coanc_admix, coanc_to_kinship')
        coancestry <- coanc_admix(admix_proportions, inbr_subpops) # the coancestry matrix
        kinship <- coanc_to_kinship(coancestry) # kinship matrix
        
        if (generations > 1) {
          # simulate realistic generations of families
          
          if (verbose)
            message('sim_children_generations_kinship')
          # defines parents semi-randomly based on avoidance of close relatives
          # but otherwise with strong assortative mating to preserve population structure
          data_G <- sim_children_generations_kinship(kinship, generations, verbose = verbose)
          parents <- data_G$parents # list with a matrix per generation
          kinship <- data_G$kinship # final kinship matrix
          
          if (verbose)
            message('sim_children_generations_admix_proportions')
          # get correct admixture proportions for these children!
          admix_proportions <- sim_children_generations_admix_proportions(admix_proportions, parents, verbose = verbose)
          
          if (verbose)
            message('sim_children_generations_genotypes')
          # final genotypes
          X <- sim_children_generations_genotypes(X, parents, verbose = verbose)
          
          # NOTE: p_anc doesn't change, this is accurate
          
          # handle fixed loci (a big pain!)
          fixed_loci_indexes <- fixed_loci(X)
          m_loci_fixed <- sum( fixed_loci_indexes )
          while (m_loci_fixed > 0) { # draw things anew, over and over until nothing was fixed
            # draw allele freqs and genotypes
            out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci_fixed)
            # overwrite fixed loci with redrawn polymorphic data
            #            X[fixed_loci_indexes, ] <- out$X # genotypes
            p_anc[fixed_loci_indexes] <- out$p_anc # ancestral AFs
            
            if (verbose)
              message('sim_children_generations_genotypes (redrawn)')
            # repeat children draws through generations, then
            # overwrite fixed loci with redrawn (hopefully) polymorphic data
            X[fixed_loci_indexes, ] <- sim_children_generations_genotypes(out$X, parents, verbose = verbose)
            
            # look for remaining fixed loci (to continue or stop loop)
            fixed_loci_indexes <- fixed_loci(X)
            m_loci_fixed <- sum( fixed_loci_indexes )
          }
        }
        
        ################
        ### 2: TRAIT ###
        ################
        
        if (verbose) message('simtrait')
        
        # parameters of simulation
        
        # create simulated trait and associated data
        # use version for known p_anc (prefered, only applicable to simulated data)
        obj_trait <- sim_trait(X = X, m_causal = m_causal, herit = herit, p_anc = p_anc)
        
        #################
        ### 3: RETURN ###
        #################
        
        # return a few items of interest
        # NOTES:
        # - we never use coancestry outside, only kinship
        # - inbr_subpops, p_anc also not used outside
        # - could return cousins_tib, don't right now (has to be conditional on cousins==TRUE)
        list(
          X = X, # genotype matrix
          kinship = kinship, # kinship matrix
          admix_proportions = admix_proportions, # admixture proportions
          trait = obj_trait$trait, # trait vector
          causal_indexes = obj_trait$causal_indexes, # randomly-picked causal locus index
          causal_coeffs = obj_trait$causal_coeffs # locus effect size vector (for causal loci only, rest are zero)
        )
      }
      
  # simulate genotypes and trait as usual
  obj <- sim_geno_trait_k3(
    n_ind = n_ind,
    m_loci = m_loci,
    m_causal = m_causal,
    k_subpops = 10,
    bias_coeff = bias_coeff,
    generations = generations,
    herit = herit,
    verbose = FALSE,
    fst = fst
  )
  X <- obj$X
  admix_proportions <- obj$admix_proportions
  kinship <- obj$kinship
  trait <- obj$trait
  causal_indexes <- obj$causal_indexes
  
  ## how we'd get these values if we weren't setting them on command line
  n_ind <- ncol(X)
  m_loci <- nrow(X)
  k_subpops <- ncol(admix_proportions)
  ## sanity checks instead
  stopifnot( ncol(X) == n_ind )
  stopifnot( nrow(X) == m_loci )
  stopifnot( ncol(admix_proportions) == k_subpops )
  
 

    # compute true eigenvectors (for PCA oracle versions)
    eigenvectors <- kinship_to_evd(kinship)

    kinship_tibble <- as_tibble( 2 * kinship, .name_repair = 'minimal')

  ########################
  ### ASSOCIATION TEST ###
  ########################
  
  # store runtimes in these tibbles
  times_kin <- tibble(.rows = 1)
  times_gas <- tibble(.rows = 1)
  pvals <- tibble(.rows = m_loci)
  
  
  ###########
  ### PCA ###
  ###########
  
  name <- "PCA"
  message(name)
  
  times_kin[[name]] <- system.time({
    kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
    eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues
  })[3]
  
  indexes <- 1 : pcs
  obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, indexes ])
  times_gas[[name]] <- obj$runtime
  pvals[[name]] <- obj$pvals
  
    
  style <- read_tsv('style.txt', col_types = 'cci')

  
  

    
  pvals_auc_collect  <- function (name_out, causal_indexes, pval_tab, times_kin, times_gas, style, cex_leg = 0.7) {
      # we shall generate other metrics in the plots below
      pval_rmsd <- list() # p-value mean distance from desired (under null)
      pr_auc <- list() # AUCs for P-R curves
      
      # reorder input tables as desired (style table)
      # also performs sanity checks (checks if style table is missing names)
      pval_tab <- reorder_tibble_names(style$name, pval_tab)
      times_kin <- reorder_tibble_names(style$name, times_kin)
      times_gas <- reorder_tibble_names(style$name, times_gas)
 
      # names and style order for p-value figs
      myNames <- names(pval_tab)
      order <- match(myNames, style$name)
      
      for (myName in myNames) {
        # sort p-values, compare to uniform quantiles, calculate RMSD
        obj <- pvals_to_null_rmsd(pval_tab[[myName]], causal_indexes)
        # store RMSD value
        pval_rmsd[[myName]] <- obj$rmsd
        # find row in style, to get col/lty from
        pr <- pvals_to_pr(pval_tab[[myName]], causal_indexes)
        # store AUCs
        pr_auc[[myName]] <- pr$auc.integral
        
     
      }
      return(c( pval_rmsd, pr_auc))
      }
  
  
  #store mean of p and auc values after removing the casual index
  p<-pvals_auc_collect(name_out, causal_indexes, pvals, times_kin, times_gas, style=style)
  p<-unlist(p)
  M_auc[j,pcs]<-p[2]
  M_rmsd[j,pcs]<-p[1]
    }
 
    
  }
  setwd("/dscrhome/yy222/fix_k/k_10")
  
  write.table(M_auc, file="auc_k_fixed_k_10_pcs_105_295_n_100.txt")
  write.table(M_rmsd, file="rmsd_k_fixed_k_10_pcs_105_295_n_100.txt")
  
  
