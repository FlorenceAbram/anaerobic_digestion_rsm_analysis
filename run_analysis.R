set.seed(20260717)

source(file.path("R", "01_reconstruct_rsm_models.R"))
source(file.path("R", "02_anova_lack_of_fit.R"))
source(file.path("R", "03_model_diagnostics.R"))
source(file.path("R", "04_coliform_and_design_space.R"))
source(file.path("R", "05_validation_uncertainty_importance.R"))

phase2a_reconstruct_main()
phase2a_anova_main()
phase2a_diagnostics_main()
phase2a1_main()
phase2b_main()
