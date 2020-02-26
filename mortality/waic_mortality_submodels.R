
#### Aims of prog: Calculate the waic for the submodels (Gelman2013, p. 181-183)

#### Load packages
library(doParallel)
library(rstanarm)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Common variables
nbFormulae = 7

#### Create the cluster
## Cluster variables
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("array id = ", array_id))

(formula_id = ifelse(array_id %% nbFormulae == 0, nbFormulae, array_id %% nbFormulae))
(species_id = array_id %/% nbFormulae + ifelse(array_id %% nbFormulae == 0, 0, 1))

nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

print("nodeslist")
print(nodeslist)
print("end nodeslist")

## Make cluster
cl = makeCluster(nodeslist, type = "PSOCK")
registerDoParallel(cl)

print("cluster done")

#### Tool functions taken from rstanarm package
source("./loo.R")
source("./misc.R")
source("./log_lik.R")
source("./stan_glm.fit.R")
source("./lArgsStanReg_submodels.R")
source("./as.matrix.stanreg.R")

my_waic <- function(x, species_id, formula_id, rm_chains = NULL, ...)
{
	if (!used.sampling(x))
		STOP_sampling_only("waic")
	if (is.stanjm(x)) {
		out <- waic.matrix(log_lik(x))
	} else if (is.stanmvreg(x)) {
		M <- get_M(x)
		ll <- do.call("cbind", lapply(1:M, function(m) log_lik(x, m = m)))
		out <- waic.matrix(ll)
	} else if (is_clogit(x)) {
		out <- waic.matrix(log_lik(x))
	} else {
		if (file.exists(paste0("array_", species_id, "/out", formula_id, "-submodels.rds")))
		{
			print(paste0("array_", species_id, "/out", formula_id, ".rds exists. Loading ..."))
			args <- readRDS(paste0("array_", species_id, "/out", formula_id, "-submodels.rds"))
		} else {
			args <- my_ll_args.stanreg(x, species_id = species_id, formula_id = formula_id, rm_chains = rm_chains)
		}
		print(paste0("Starting loo::waic for species ", species_id, " formula ", formula_id))
		out <- loo::waic.function(ll_fun(x), data = args$data, draws = args$draws)
	}
	structure(out,
		class = c("waic", "loo"),
        model_name = deparse(substitute(x)),
        discrete = is_discrete(x),
        yhash = hash_y(x),
        formula = loo_model_formula(x))
}

my_as.matrix.stanreg <- function(x, ..., rm_chains = NULL, pars = NULL, regex_pars = NULL) {
	pars <- collect_pars(x, pars, regex_pars)
	user_pars <- !is.null(pars)

	if (used.optimizing(x))
	{
	  	mat <- x$asymptotic_sampling_dist
	    if (is.null(mat))
	      STOP_no_draws()
	    if (!user_pars)
		{
	    	aux <- c("sigma", "scale", "shape", "lambda", "reciprocal_dispersion")
			pars <- c(names(coef(x)), # return with coefficients first
				aux[which(aux %in% colnames(mat))])
    	}
	} else { # used mcmc or vb
		mat <- as.matrix(x$stanfit) # I will remove the non-converged chains of this object
		## Remove the chains I do not want to keep:
		# The matrix mat is organised as folow:
		#	- nb_iterToConsider = nb_iter after thinning
		#	- nrow = nb_chains * nb_iterToConsider
		#	- each column is a parameters
		#	- it is recursively organised first chain = 1:nb_iterToConsider
		#		second chain = (nb_iterToConsider + 1):(2*nb_iterToConsider)
		#		k^th chain = ((k-1)nb_iterToConsider + 1):(k*nb_iterToConsider)

		if (!is.null(rm_chains))
		{
			# Check the number of chains & iterations, for simplicity, use an even number of iterations
			nb_chains = x$stanfit@sim[["chains"]]
			nb_iterToConsider = x$stanfit@sim[["iter"]] - x$stanfit@sim[["warmup"]]
			if (min(rm_chains) < 0 | max(rm_chains) > nb_chains)
				print(paste0("*** Error, rm_chains must be positive, and lower than ", nb_chains))
			rm_chains = unique(rm_chains)
			nb_to_del = nb_iterToConsider*ncol(mat)
			for (i in 1:length(rm_chains))
			{
				kth = rm_chains[i] # the k^th chain
				mat[((kth - 1)*nb_iterToConsider + 1):(kth*nb_iterToConsider), ] = rep(NA, nb_to_del)
			}
			mat = matrix(data = mat[!is.na(mat)], nrow = nrow(mat) - nb_iterToConsider*length(rm_chains))
		}

	    if (!user_pars)
	    	pars <- exclude_lp_and_ppd(colnames(mat))
	}
	if (user_pars)
		check_missing_pars(mat, pars)

	# ncol = nb_explanatory + intercept + nb_randEff + Sigma[plot_id:(Intercept),(Intercept)]
	# 	 = 20 + 1 + nbPlot + 1 in my case
	mat <- mat[, pars, drop = FALSE]
	if (!is.mer(x))
		return(mat)
	unpad_reTrms(mat) # It removes "b[(Intercept) plot_id:_NEW_plot_id]", which I do not know what it is
}

get_z.lmerMod <- function(object, ...)
{
	Zt <- object$glmod$reTrms$Zt %ORifNULL% stop("Z not found")
	Matrix::t(Zt)
}

#### Load data and run WAIC
## Model
(loadPath = paste0("./array_", species_id, "/"))
(ls_files = list.files(path = loadPath, pattern = "submodel_"))

model = readRDS(paste0(loadPath, ls_files[formula_id]))

waic_model = my_waic(x = model, species_id = species_id, formula_id = formula_id, rm_chains = NULL)
saveRDS(waic_model, paste0(loadPath, "waic", formula_id, "_submodels.rds"))
