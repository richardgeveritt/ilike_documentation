# Compiles every example of the book's language chapters and runs each one
# briefly, so that the book cannot drift from the package. Run from this
# folder: Rscript check_examples.R
library(ilike)
out = tempfile("ilike_examples_")
dir.create(out)
run = list(
  first.ilike = function(r) MCMC(r, results_path = out, results_name = "first", initial_values = list(list(theta = c(0, 0))),
                                 model_parameter_list = list(proposal_covariance = 0.01 * diag(2)), seed = 1),
  schools.ilike = function(r) MCMC(r, results_path = out, results_name = "schools", seed = 1),
  regression.ilike = function(r) MCMC(r, results_path = out, results_name = "regression", seed = 1),
  ising.ilike = function(r) MCMC(r, results_path = out, results_name = "ising", seed = 1),
  smooth.ilike = function(r) MCMC(r, results_path = out, results_name = "smooth", seed = 1),
  user_functions.ilike = function(r) MCMC(r, algorithm = "mcmc", results_path = out, results_name = "user", seed = 1),
  r_functions.ilike = function(r) MCMC(r, results_path = out, results_name = "r", seed = 1),
  abc.ilike = function(r) { SMC_with_MCMC(r, algorithm = "abc_smc", results_path = out, results_name = "abc_smc", seed = 1)
                            MCMC(r, algorithm = "mcmc", results_path = out, results_name = "abc", seed = 1,
                               initial_values = list(list(theta = 1))) },
  pmmh.ilike = function(r) MCMC(r, results_path = out, results_name = "pmmh", seed = 1,
                                initial_values = list(list(a = 0.5))),
  enki.ilike = function(r) IENKI(r, results_path = out, results_name = "enki", seed = 1),
  ssm.ilike = function(r) { KF(r, algorithm = "kf", results_path = out, results_name = "kf")
                            PF(r, algorithm = "pf", results_path = out, results_name = "pf", seed = 1)
                            EnKF(r, algorithm = "enkf", results_path = out, results_name = "enkf", seed = 1) },
  ssm_guided.ilike = function(r) { PF(r, algorithm = "bootstrap", results_path = out, results_name = "bootstrap", seed = 1)
                                   PF(r, algorithm = "guided", results_path = out, results_name = "guided", seed = 1) },
  mfds.ilike = function(r) MFDS(r, results_path = out, results_name = "mfds", seed = 1),
  ssm_inner.ilike = NULL  # (compiled only: its data and a come from pmmh.ilike)
)
files = list.files(".", pattern = "\\.ilike$")
missing = setdiff(files, c(names(run), "ssm_inner.ilike"))
if (length(missing) > 0) stop("no run given for: ", paste(missing, collapse = ", "))
for (f in files)
{
  cat(f, ": ", sep = "")
  recipe = compile(f)
  if (!is.null(run[[f]])) run[[f]](recipe)
  cat("ok\n")
}
unlink(out, recursive = TRUE)
