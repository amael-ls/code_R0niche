my_ll_args.stanreg <- function(object, newdata = NULL, offset = NULL, m = NULL,
	reloo_or_kfold = FALSE, species_id, formula_id, rm_chains = NULL, ...)
{
	print("my_ll_args.stanreg starting")
	validate_stanreg_object(object)
	f <- family(object, m = m)
	draws <- nlist(f)
	has_newdata <- !is.null(newdata)

	dots <- list(...)

	z_betareg <- NULL
	stanmat <- my_as.matrix.stanreg(object, rm_chains)
	x <- get_x(object, m = m)
	y <- get_y(object, m = m)

	if (!is_polr(object))
	{
		fname <- f$family

	    trials <- 1
	    if (is.factor(y))
	    	y <- fac2bin(y)
	    stopifnot(all(y %in% c(0, 1)))
		    data <- data.frame(y, trials, x)
    }
    nms <- if (is.stanmvreg(object))
    	collect_nms(colnames(stanmat),
	  		M = get_M(object),
			stub = get_stub(object)) else NULL
    beta_sel <- if (is.null(nms)) seq_len(ncol(x)) else nms$y[[m]]
    draws$beta <- stanmat[, beta_sel, drop = FALSE]
    m_stub <- get_m_stub(m, stub = get_stub(object))

	data$offset <- if (has_newdata) offset else object$offset
	if (is.mer(object))
	{
		b_sel <- if (is.null(nms)) b_names(colnames(stanmat)) else nms$y_b[[m]]
		b <- stanmat[, b_sel, drop = FALSE]
		z <- get_z(object, m = m)
		if (unique(z@x != 1))
			print("*** ERROR from my_ll_args.stanreg ***: I assumed the matrix to have only 0s and 1s")

		if (unique(z@x == 1))
		{
			print("There are only 0 and 1 in the dense matrix")
			Matrix::writeMM(z, paste0("array_", species_id, "/z", formula_id, "-submodels.txt"))
			# Position of the 1s (i.e., the only non-zero coeff of the matrix)
			pos_1 = data.table::fread(paste0("array_", species_id, "/z", formula_id, "-submodels.txt"), header = FALSE)
			n_row = nrow(z)
			n_col = ncol(z)
			z2 = matrix(data = 0, nrow = n_row, ncol = n_col)
			z2_ind = (pos_1[, V2] - 1)*n_row + pos_1[, V1]
			z2[z2_ind] = 1
			test = z2 - z
			if (unique(test@x) == 0)
				print("transformation succeded")
		}
    }
    data <- cbind(data, z2[1:NROW(x),, drop = FALSE])
    draws$beta <- cbind(draws$beta, b)
    out <- nlist(data, draws, S = NROW(draws$beta), N = nrow(data))
	saveRDS(out, paste0("array_", species_id, "/out", formula_id, "-submodels.rds"))
	print("lArg done")
	return(out)
}
