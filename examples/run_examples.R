# Runs the examples of the "Examples" chapter and saves the output that the
# chapter plots (in output/). Run from this directory:
#   Rscript run_examples.R [section ...]
# where section is one or more of gaussian, smc, cv, lv_mcmc, lv_abc_smc
# (default: all). The chapter does not run ilike itself.

library(ilike)

sections = commandArgs(trailingOnly = TRUE)
if (length(sections) == 0)
  sections = c("gaussian", "smc", "cv", "lv_mcmc", "lv_abc_smc")

results = file.path(tempdir(), "ilike_examples")
dir.create(results, showWarnings = FALSE)

# keep every k-th iteration of each chain, to keep the stored output small
thin = function(output, k) output[output$Iteration %% k == 0, ]

if ("gaussian" %in% sections)
{
  timings = c()
  for (algorithm in c("rwm", "langevin", "barker"))
  {
    recipe = compile("gaussian_model.ilike")
    time = system.time(MCMC(recipe, algorithm = algorithm, results_name = algorithm, results_path = results,
                            initial_values = list(list(θ = 2)), seed = 1))[["elapsed"]]
    output = load_mcmc_output(file.path(results, algorithm))
    output$Algorithm = algorithm
    saveRDS(output, file.path("output", paste("gaussian_", algorithm, ".rds", sep = "")))
    timings[algorithm] = time
  }
  saveRDS(timings, file.path("output", "gaussian_timings.rds"))
}

if ("smc" %in% sections)
{
  recipe = compile("gaussian_model.ilike")
  log_evidence = SMC_with_MCMC(recipe, algorithm = "smc", results_name = "smc", results_path = results,
                               number_of_particles = 1000, smc_iterations_to_store = 1000, seed = 1)
  output = load_smc_output(file.path(results, "smc"))
  # keep the last MCMC iteration of each particle at each target
  output = output[output$Iteration == ave(output$Iteration, output$Target, FUN = max), ]
  saveRDS(list(log_evidence = log_evidence, output = output), file.path("output", "gaussian_smc.rds"))
}

if ("cv" %in% sections)
{
  pl = list(prior_mean = c(0, 0), prior_covariance = 100*diag(2))
  model = compile("cv_model.ilike")

  kf = KF(model, algorithm = "kf", results_name = "cv_kf", results_path = results, model_parameter_list = pl,
          kf_iterations_to_store = 21)
  pf_estimates = function(algorithm)
    sapply(1:50, function(s) PF(model, algorithm = algorithm, results_name = paste("cv", algorithm, s, sep = "_"),
                                 results_path = results, model_parameter_list = pl, number_of_particles = 100,
                                 seed = s))
  bootstrap = pf_estimates("bootstrap")
  guided_estimates = pf_estimates("guided")

  pf = PF(model, algorithm = "guided", results_name = "cv_pf", results_path = results, model_parameter_list = pl,
          number_of_particles = 1000, smc_iterations_to_store = 21, seed = 1)
  output = load_smc_output(file.path(results, "cv_pf"))
  PF(model, algorithm = "bootstrap", results_name = "cv_pf_bootstrap", results_path = results,
     model_parameter_list = pl, number_of_particles = 1000, smc_iterations_to_store = 21, seed = 1)
  bootstrap_output = load_smc_output(file.path(results, "cv_pf_bootstrap"))
  ess = rbind(data.frame(Filter = "Bootstrap", unique(bootstrap_output[, c("Target", "ISESS")])),
              data.frame(Filter = "Guided", unique(output[, c("Target", "ISESS")])))

  utils::data("constant_velocity_x", package = "ilike")
  saveRDS(list(kf = kf, bootstrap = bootstrap, guided = guided_estimates, pf_output = output, ess = ess,
               truth = constant_velocity_x[1:20, ]),
          file.path("output", "cv.rds"))
}

if ("lv_mcmc" %in% sections)
{
  lv = list()
  for (epsilon in c(1, 0.5))
  {
    pl = list(number_of_simulations = 10, tolerance = epsilon, proposal_covariance = 0.01*diag(3))
    recipe = compile("lv_model.ilike")
    name = paste("lv_abc_mcmc_", epsilon, sep = "")
    time = system.time(MCMC(recipe, algorithm = "abc_mcmc", results_name = name, results_path = results, model_parameter_list = pl,
                            initial_values = list(list(log_theta = log(c(1, 0.005, 0.6)))), seed = 1))[["elapsed"]]
    output = thin(load_mcmc_output(file.path(results, name)), 10)
    output$Epsilon = epsilon
    lv[[as.character(epsilon)]] = list(output = output, time = time)
  }
  saveRDS(lv, file.path("output", "lv_abc_mcmc.rds"))
}

if ("lv_abc_smc" %in% sections)
{
  # (the tolerance constant is not used: the SMC sets epsilon)
  pl = list(number_of_simulations = 10, proposal_covariance = 0.01*diag(3))
  recipe = compile("lv_model.ilike")
  time = system.time(log_evidence <- SMC_with_MCMC(recipe, algorithm = "abc_smc", results_name = "lv_abc_smc", results_path = results,
                                                   model_parameter_list = pl, number_of_particles = 1000,
                                                   smc_iterations_to_store = 1000, seed = 1))[["elapsed"]]
  output = load_smc_output(file.path(results, "lv_abc_smc"))
  # keep the state of each particle after its last MCMC iteration at each target
  output = output[output$Iteration == ave(output$Iteration, output$Target, FUN = max), ]
  saveRDS(list(log_evidence = log_evidence, output = output, time = time), file.path("output", "lv_abc_smc.rds"))
}
