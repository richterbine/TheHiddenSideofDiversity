# TheHiddenSideofDiversity

Studies on ecological communities often address patterns of species distribution and abundance, but few consider uncertainty in counts of both species and individuals when computing diversity measures. We evaluated the extent to which imperfect detection may influence patterns of taxonomic, functional and phylogenetic diversity in ecological communities. For this, we employed a hierarchical Bayesian N-mixture model to estimate the true abundance of fruit-feeding butterflies sampled in canopy and understory strata in a subtropical forest. We then developed a tool to calculate how much information is hidden when uncertainty in counts of species and individuals is not considered when estimating diversity.

Descriptions for the files and scripts

The `R` folder contains.

1. Bayesian_models
    1. *Static_Nmixture_P_RE.R* - script to generate the Bayesian model
    1. *Static_Nmixture_P_RE.txt* - the model code in a .txt file.

1. Model_validation
   1. *D_Model_validation.R* -  script to generate the simulated communities and perform the bayesian model.
   1. *S_Model_sim_evaluation.R* - script to evaluate ability of the bayesian model in found the real parameters.
   1. *Comm_simulated.rds* - the simulated communities used to run the bayesian model.
   1. *model_comm_simulated.rds* - output of the model for each simulated community. 

1. functions 
    1. *hidden_diversity_function.R* - function to calculate the hidden diversity for each diversity measures. 
    2. *hidden_diversity_parallel.R* - function parallelized, a faster version of hidden diversity (in construction).
    3. *simCommNmix_function.R* - fucntion to generate communities, accounting with the imperfect detection.
    4. *genus_to_spp_function.R* - function to manipulated and prune a phylogenetic tree to a specific species pool.

1. *M_Data_manipulations.R* - script to read the community data observed, the phylogenetic tree and the traits matrix.
1. *D_S_Model_bfly_data.R* - script to run the bayesian model with the empirical dataset (fruit-feeding butterflies community sampled at FLONA-SFP between 2016 and 2017).
1. *D_Hidden_diversity_bfly_data.R* - script to perform the hidden diversity for the empirical dataset.
2. *S_HD_bfly.R* - script to plot and analyse the hidden diversity framework.

The `output` folder contains:

1. figures
    1. *Fig1_HD_framework.png* - Figure 1 main text.
    1. *Fig2_meanC.png* - Figure 2 main text.
    2. *Fig3_HD.png* - Figure 3 main text.
    3. *FigS1_meanCparams.png* - Supplementary material Figure S1.
    4. *FigS2_sdCparams.png* - Supplementary material Figure S2.
    5. *FigS3_A-F_hyperparams.png* - Supplementary material Figure S3.
    6. *FigS4_G-L_hyperparams.png* - Supplementary material Figure S4.
    7. *FigS5_meanSparams.png* - Supplementary material Figure S5.
    8. *FigS6_pattern_obs-est.png* - Supplementary material Figure S6.
    
1. *out.mod.bfly.rds* - the bayesian model output with empirical data.
1. *fitjags.mod.bfly.rds* - output for the basic.jags model, with the MCMC for true abundance matrix.
1. *joint_bfly_data.txt* - data frame with merged data for temperature and species counts of the empirical dataset.
2. *HD_bfly.rds* - output with 100 estimated values for each diversity measure + diversity values calculated with observed data.
3. *tree_bfly_FLONA.txt* - the pruned phylogenetic tree for the empirical dataset.
4. *tree.func_bfly_FLONA.txt* - the functional dendrogram based on similarity of functional traits among species.
5. *mod_lmer.txt* - output with intercepts and slopes for linear mixed models between HD and strata (canopy x undesrtory).

The `data` folder contains:
1. processed
    1. *Mean_traits_bfly.txt* - the functional traits matrix of the empirical dataset.
    1. *raw_data_bfly.txt* - the raw count data of t.
    1. *raw_temp_bfly.txt* - the measures of the temperature for each trap, in each site by sampling days in each strata.

1. Raw
    1. *Nymphalidae_tree.txt* - the phylogenetic tree proposed by Chazot et al. 2019.
    1. *Raw_data_FLONA.xlsx* - the raw data for fruit-feeding butterflies communities sampled at FLONA-SFP.
