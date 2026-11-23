#' @title
#' Simulate continuous data from a multilevel graphical VAR model
#'
#' @description
#' Simulates continuous multivariate time series for one or more persons
#' from a multilevel graphical VAR(1) model. The simulation uses several
#' functions adapted from the `modnets` package to generate data from
#' specified or generated temporal, contemporaneous, and between-person
#' network structures.
#'
#' @param ntime Number of retained observations per person. A scalar is repeated.
#' @param nperson Number of persons.
#' @param nnode Number of nodes.
#' @param lag Number of lags. This reduced version supports `lag = 1`.
#' @param temporal_sparsity Proportion of all `nnode^2` population temporal
#'   coefficients set exactly to zero, including autoregressive coefficients.
#' @param temporal_prop_positive Probability that a nonzero population temporal
#'   coefficient is positive.
#' @param beta_fixed Optional stationary `nnode` x `nnode` matrix of
#'   population temporal VAR(1) regression coefficients. `NULL` generates it.
#'   When supplied, `temporal_sparsity` and `temporal_prop_positive` are ignored
#'   with a message; temporal random-effect settings remain in effect.
#' @param temporal_beta_re_sd_range Lower and upper bounds for the random-effect
#'   SDs of person-specific temporal regression coefficients; default `c(.1, 1)`.
#' @param temporal_re_correlation_sparsity Proportion of zero off-diagonal
#'   marginal correlations among temporal random slopes. A value of `1`
#'   makes these slopes mutually uncorrelated. Mean-slope correlations are zero.
#' @param contemporaneous One of `"wishart"`, `"random_ggm"`, `"orthogonal_re"`
#'   or `"fixed"`. For `"orthogonal_re"`, independent deviations are applied
#'   to unique contemporaneous partial-correlation edges while the precision
#'   diagonal is fixed.
#' @param contemporaneous_sparsity Proportion of zero edges in the generated
#'   contemporaneous partial-correlation network. For `contemporaneous = "random_ggm"`,
#'   it applies to each person-specific network rather than to
#'   the population network. Ignored if `theta_fixed` is supplied.
#' @param theta_var Positive vector of `nnode` innovation variances; default
#'   `rep(1, nnode)`. If `theta_fixed` is provided, `theta_var` is ignored and
#'   replaced by `diag(theta_fixed)` with a message
#' @param theta_fixed Optional positive-definite `nnode` x `nnode`
#'   innovation covariance matrix. Overrides `contemporaneous_sparsity` and
#'   `theta_var` with a message Supplying it for `"random_ggm"` causes an
#'   error because its population covariance is computed from individual draws.
#' @param contemporaneous_pcor_re_sd_range Lower and upper bounds for the
#'   random-effect SDs of unique contemporaneous partial-correlation edges.
#'   Used only for `contemporaneous = "orthogonal_re"`.
#' @param between_sparsity Sparsity of the between-person partial-correlation
#'   network formed by person-specific node means.
#' @param fixed_mean_sd SD used to draw population node means from
#'   `N(0, fixed_mean_sd^2)`.
#' @param random_intercept_sd_range Lower and upper bounds for the random-effect
#'   SDs of person-specific node means; default `c(1, 1)`.
#' @param between_cov Optional positive-definite `nnode` x `nnode`
#'   covariance matrix of the random intercepts. If supplied, its diagonal
#'   fixes the random-intercept variances; `between_sparsity` and
#'   `random_intercept_sd_range` are ignored with a message.
#' @param shrink_temporal_fixed Whether population temporal coefficients may
#'   shrink following a nonstationary individual network.
#' @param shrink_temporal_fixed_factor Multiplicative factor for shrinking
#'   population temporal coefficients when enabled.
#' @param shrink_temporal_re_factor Multiplicative factor for shrinking
#'   temporal random-slope SDs after failed temporal draws.
#' @param shrink_contemporaneous_re_factor Multiplicative factor for shrinking
#'   contemporaneous partial-correlation random-effect SDs if an individual
#'   precision matrix fails positive-definiteness checks.
#' @param max_shrink_attempts Maximum attempts for each shrinking loop.
#' @param only_nets If `TRUE`, return population and person-specific model
#'   matrices without simulating the time-series observations.
#'
#' @return An `mlGVARsim`-like object. `model_parameter$generation_info` reports requested and
#'   realized sparsities and the final random-effect SDs after possible shrinking.
#'
#' @details
#'
#' Compared with `modnets::mlGVARsim()`, this function omits moderation
#' effects. It retains person-specific means, temporal coefficients, and
#' contemporaneous innovation covariance matrices, as well as a burn-in
#' period. When `beta_fixed` is generated, candidate population matrices
#' are drawn until the VAR(1) stationarity condition is met.
#'
#' The function first sets up the population parameters. If `beta_fixed`
#' is generated, candidate temporal coefficient matrices are drawn until
#' a stationary matrix is obtained. Person-specific innovation covariances
#' are determined by `contemporaneous`. For `"random_ggm"`, the population
#' innovation covariance is defined as the average of the generated
#' person-specific covariance matrices.
#'
#' The function then draws person-specific means and temporal coefficients.
#' If any person's temporal matrix is nonstationary, the function reduces
#' the temporal random-slope SDs and redraws the person-specific parameters.
#' If `shrink_temporal_fixed = TRUE`, it also scales down the population
#' temporal coefficients. For `"orthogonal_re"`, a separate step reduces
#' the contemporaneous random-effect SDs if a person's precision matrix
#' fails the positive-definiteness check.
#'
#' Each time series starts at the person's node means. The first 100
#' observations are discarded as burn-in, leaving `ntime` observations
#' per person. With `only_nets = TRUE`, the function returns the generating
#' parameters and networks before simulating observations.
#'
#' @export
mlvar_simulate <- function(
  # Study design
  ntime = 50, nperson = 10, nnode = 3, lag = 1,
  # Temporal population network and its random effects
  temporal_sparsity = 0.5,
  temporal_prop_positive = 0.5,
  beta_fixed = NULL,
  temporal_beta_re_sd_range = c(0.1, 0.5),
  temporal_re_correlation_sparsity = 0.5,
  # Contemporaneous network and its random effects
  contemporaneous = c("wishart", "random_ggm", "fixed", "orthogonal_re"),
  contemporaneous_sparsity = 0.5,
  theta_var = rep(1, nnode),
  theta_fixed = NULL,
  contemporaneous_pcor_re_sd_range = c(0.05, 0.05),
  # Between-person network and random intercepts
  between_sparsity = 0.5,
  fixed_mean_sd = 1,
  random_intercept_sd_range = c(1, 1),
  between_cov = NULL,
  # Shrinking controls
  shrink_temporal_re_factor = 0.9,
  shrink_contemporaneous_re_factor = 0.9,
  shrink_temporal_fixed = FALSE,
  shrink_temporal_fixed_factor = 0.9,
  max_shrink_attempts = 1000L,
  # Output mode
  only_nets = FALSE) {

  # Select the contemporaneous method; no legacy names are supported.
  contemporaneous <- match.arg(contemporaneous)

  if (!is.numeric(lag) || length(lag) != 1L || !isTRUE(lag == 1)) {
    stop("mlvar_simulate() currently supports lag = 1 only.", call. = FALSE)
  }

  # Validate arguments (including defaults); the validator never modifies them.
  # Simulation dimensions and output
  checker(
    nperson = nperson,
    ntime = ntime,
    nnode = nnode,
    only_nets = only_nets
  )

  # Temporal network
  checker(
    nnode = nnode,
    beta_fixed = beta_fixed,
    temporal_sparsity = temporal_sparsity,
    temporal_prop_positive = temporal_prop_positive,
    temporal_beta_re_sd_range = temporal_beta_re_sd_range,
    temporal_re_correlation_sparsity = temporal_re_correlation_sparsity
  )

  # Contemporaneous network
  checker(
    nnode = nnode,
    contemporaneous = contemporaneous,
    theta_fixed = theta_fixed,
    theta_var = theta_var,
    contemporaneous_sparsity = contemporaneous_sparsity,
    contemporaneous_pcor_re_sd_range = contemporaneous_pcor_re_sd_range
  )

  # Between-person network and means
  checker(
    nnode = nnode,
    between_cov = between_cov,
    between_sparsity = between_sparsity,
    random_intercept_sd_range = random_intercept_sd_range,
    fixed_mean_sd = fixed_mean_sd
  )

  # Shrinkage
  checker(
    shrink_temporal_fixed = shrink_temporal_fixed,
    shrink_temporal_fixed_factor = shrink_temporal_fixed_factor,
    shrink_temporal_re_factor = shrink_temporal_re_factor,
    shrink_contemporaneous_re_factor = shrink_contemporaneous_re_factor,
    max_shrink_attempts = max_shrink_attempts
  )

  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop(
      "Package 'mvtnorm' is required to use this function.",
      call. = FALSE
    )
  }


  # Expand a shared observation count; parameter defaults are set in the signature.
  if (length(ntime) == 1L) ntime <- rep(ntime, nperson)



  # Explicit population matrices take priority over their generation settings.
  # Warn at the public entry point; the validator above performs error checks only.
  if (!is.null(beta_fixed)) {
    message("beta_fixed supplied: temporal_sparsity and ",
            "temporal_prop_positive are ignored.")
  }
  if (!is.null(theta_fixed)) {
    message("theta_fixed supplied: contemporaneous_sparsity and theta_var ",
            "are ignored; theta_var is set to diag(theta_fixed).")
    theta_var <- diag(theta_fixed)
  }
  if (!is.null(between_cov)) {
    message("between_cov supplied: between_sparsity and ",
            "random_intercept_sd_range are ignored.")
  }
  # The edge-level partial-correlation SDs are used only by orthogonal_re.
  if (contemporaneous != "orthogonal_re" &&
      !missing(contemporaneous_pcor_re_sd_range)) {
    warning("contemporaneous_pcor_re_sd_range is ignored when contemporaneous = '",
            contemporaneous, "'.", call. = FALSE)
  }

  # Model dimensions derived from nnode.
  n_temporal <- nnode^2
  n_random_effects <- nnode + n_temporal

  # Record which matrices were explicitly supplied to distinguish given from generated sparsity.
  beta_fixed_input <- beta_fixed
  between_cov_input <- between_cov
  theta_fixed_input <- theta_fixed

  # Define all independently specified population parameters first:
  # random-effect distribution, contemporaneous network, mu_fixed and
  # stationary beta_fixed. Exception: under random_ggm, theta_fixed is the
  # average of the subsequently generated person-specific covariances.
  # Draws and their order follow the previous refactored implementation.
  population <- setup_population(
    nnode = nnode,
    n_temporal = n_temporal,
    n_random_effects = n_random_effects,
    random_intercept_sd_range = random_intercept_sd_range,
    temporal_beta_re_sd_range = temporal_beta_re_sd_range,
    between_sparsity = between_sparsity,
    temporal_re_correlation_sparsity = temporal_re_correlation_sparsity,
    theta_var = theta_var,
    contemporaneous = contemporaneous,
    contemporaneous_sparsity = contemporaneous_sparsity,
    contemporaneous_pcor_re_sd_range = contemporaneous_pcor_re_sd_range,
    fixed_mean_sd = fixed_mean_sd,
    temporal_sparsity = temporal_sparsity,
    temporal_prop_positive = temporal_prop_positive,
    beta_fixed = beta_fixed,
    between_cov = between_cov,
    theta_fixed = theta_fixed
  )

  # Keep population parameters in one object so shrinking cannot leave
  # separate local copies with outdated values.

  # Fixed means and temporal coefficients are defined BEFORE person draws.
  # The random_ggm branch computes the population innovation covariance later.

  # 2. Draw person-specific contemporaneous innovation covariances.
  # For orthogonal_re, any non-positive-definite precision matrix rejects
  # the whole batch. Shrink ONLY the contemporaneous edge RE SDs, then
  # redraw ALL persons; the fixed covariance and precision diagonal remain fixed.
  # For wishart, random_ggm, and fixed, this loop succeeds on its first pass.
  n_contemporaneous_shrink <- 0L
  repeat {
    contemporaneous_subjects <- generate_contemporaneous_subjects(
      population = population,
      nperson = nperson,
      nnode = nnode,
      contemporaneous = contemporaneous,
      contemporaneous_sparsity = contemporaneous_sparsity,
      theta_var = theta_var
    )

    # Accept the batch only when all individual precision matrices are valid.
    if (contemporaneous_subjects$valid) break

    # Identical to the original orthogonal_re shrinking step: changing these
    # SDs also changes the population distribution of contemporaneous effects.
    n_contemporaneous_shrink <- n_contemporaneous_shrink + 1L
    if (n_contemporaneous_shrink >= max_shrink_attempts) {
      stop("Could not generate positive-definite individual contemporaneous ",
           "networks within max_shrink_attempts.")
    }
    population$contemporaneous_pcor_re_sd_current <-
      population$contemporaneous_pcor_re_sd_current *
      shrink_contemporaneous_re_factor
  }

  # random_ggm defines its population covariance as the average over the
  # generated individual covariance matrices. For other methods this is
  # identical to the already specified population$theta_fixed.
  population$theta_fixed <- contemporaneous_subjects$theta_fixed

  # Simulate individual VAR(1) series, with temporal shrinking

  # This is a SECOND, separate retry loop. It keeps the already accepted
  # contemporaneous covariances fixed across all temporal retries.
  # Count temporal failures from nonstationary individual matrices.
  n_temporal_shrink <- 0

  # Repeat person-specific mean and temporal draws; if all individual
  # temporal networks are stationary, simulate their time series.
  # A failed stationarity OR data-magnitude check applies the SAME shrink
  # step below before resampling all individual means and temporal slopes.
  repeat {

    # Draw person-specific means and temporal matrices; check all matrices.
    temporal_subjects <- generate_temporal_subjects(
      nperson = nperson,
      nnode = nnode,
      mu_fixed = population$mu_fixed,
      beta_fixed = population$beta_fixed,
      omega = population$omega
    )

    # Only draw observations if ALL person-specific temporal matrices are stationary.
    if (temporal_subjects$stationary) {
      # For one person, use population parameters rather than individual draws.
      if (nperson == 1L) {
        warning("nperson = 1: person-specific parameters equal population parameters")
        temporal_subjects$person_means <- matrix(population$mu_fixed, nrow = 1)
        temporal_subjects$person_betas <- array(
          matrix(population$beta_fixed, nnode, nnode),
          c(nnode, nnode, 1)
        )
        contemporaneous_subjects$theta <- array(population$theta_fixed, c(nnode, nnode, 1))
        if (contemporaneous == "orthogonal_re") {
          contemporaneous_subjects$gamma_theta_subject_true <-
            list(contemporaneous_subjects$gamma_theta_fixed_true)
        }
      }

      # Parameter-only mode: all population and person-specific matrices
      # have been generated and checked. Skip time-series simulation and
      # proceed directly to assembling the full parameter output below.
      if (only_nets) break

      # Generate a separate time series for each person using person-specific values.
      data_list <- lapply(seq_len(nperson), function(p) {

        #  Simulate a multivariate VAR(1) time series for person p.
        # Simulation starts from the values specified in init.
        # The first 100 simulated observations are discarded
        # as burn-in to reduce dependence on the initial values.
        res <- simulate_var1(
          beta = matrix(temporal_subjects$person_betas[, , p], nnode, nnode),
          means = temporal_subjects$person_means[p, ],
          theta = matrix(contemporaneous_subjects$theta[, , p], nnode, nnode),
          ntime = ntime[p],
          init = temporal_subjects$person_means[p, ],
          burnin = 100L
        )
        # Assign variable names to the columns of the simulated time series.
        colnames(res) <- paste0("V", seq_len(nnode))

        # Convert the simulated observations into a data frame.
        res <- as.data.frame(res)
        # Add the ID of the current person.
        res$ID <- p
        # Return the simulated data for person p. lapply() collects these data frames in data_list.
        res
      })

      # Combine all person-specific data frames by rows.
      data <- do.call(rbind, data_list)
      # Check whether any simulated network variable are infinite
      # and if they are - stop the process
      if (any(!is.finite(as.matrix(data[, seq_len(nnode), drop = FALSE])))) {
        stop("Simulation produced non-finite observations.")
      }
      # absolute value greater than 100 represent very large value with means at 0
      # however instead of just quitting as in modnets only return a value as the means
      # could be for example 100, which would make this a normal value
      if (any(abs(data[, seq_len(nnode), drop = FALSE]) > 100)) {
        warning("Simulation produced observations with absolute values above 100.
                Check creating values to make sure that these are realistic values.")
      }
      break
    }

    # Common failure handler, reached when:
    # (a) at least one person-specific temporal matrix is nonstationary, OR
    # (b) at least one simulated observation has absolute value > 100.
    # Always shrink temporal RE SDs; shrink population$beta_fixed only if
    # requested. Random-intercept SDs and contemporaneous covariances remain unchanged.
    # The next repeat iteration redraws all person means and temporal coefficients.
    n_temporal_shrink <- n_temporal_shrink + 1
    if (n_temporal_shrink >= max_shrink_attempts) {
      stop("Could not generate stationary individual VAR networks and ",
           "acceptable time series within max_shrink_attempts.")
    }
    # Fixed temporal effects can remain unchanged while temporal RE SDs shrink.
    if (shrink_temporal_fixed) {
      population$beta_fixed <- population$beta_fixed * shrink_temporal_fixed_factor
    }
    beta_positions <- nnode + seq_len(n_temporal)
    population$intercept_temporal_re_sd_current[beta_positions] <-
      population$intercept_temporal_re_sd_current[beta_positions] * shrink_temporal_re_factor
    # Rebuild joint RE covariance from the updated population SDs;
    # the supplied between-person block is preserved exactly.
    population$omega <- diag(population$intercept_temporal_re_sd_current, n_random_effects, n_random_effects) %*%
      population$omega_cor %*%
      diag(population$intercept_temporal_re_sd_current, n_random_effects, n_random_effects)
    if (!is.null(between_cov_input)) {
      population$omega[seq_len(nnode), seq_len(nnode)] <- between_cov_input
    }
  }

  # Warn when the requested temporal random-effect SDs had to be reduced.
  if (n_temporal_shrink > 0) {
    total_factor <- shrink_temporal_re_factor^n_temporal_shrink
    warning(
      "Temporal random-effect SDs were shrunk ", n_temporal_shrink,
      " times (total factor ", signif(total_factor, 3), ") to obtain ",
      "stationary individual networks. The realized SDs are smaller than ",
      "requested; see model_parameters$generation_info$random_effect_sd. ",
      "To avoid this, consider a smaller temporal_beta_re_sd_range.",
      call. = FALSE
    )
  }

  if (n_contemporaneous_shrink > 0) {
    warning(
      "Contemporaneous edge random-effect SDs were shrunk ",
      n_contemporaneous_shrink, " times (total factor ",
      signif(shrink_contemporaneous_re_factor^n_contemporaneous_shrink, 3),
      ") to obtain positive-definite individual networks. ",
      "See model_parameters$generation_info$random_effect_sd.",
      call. = FALSE
    )
  }

  # A supplied beta_fixed should only be changed if the user explicitly allowed it.
  if (shrink_temporal_fixed && n_temporal_shrink > 0 &&
      !is.null(beta_fixed_input)) {
    warning(
      "The supplied beta_fixed was shrunk to obtain stationary individual networks
      (total factor ", signif(shrink_temporal_fixed_factor^n_temporal_shrink, 3),
      ") because shrink_temporal_fixed = TRUE.",
      call. = FALSE
    )
  }

  # Output

  # Save the FINAL parameters actually used (possibly after
  # shrinking), plus the realized individual parameter draws. This allows
  # estimators to be compared against their generating values.
  model_parameters <- list(
    # Person-specific node means (random intercepts).
    mu = make_parameter_output(
      mean = population$mu_fixed,
      SD = population$intercept_temporal_re_sd_current[seq_len(nnode)],
      subject = lapply(seq_len(nperson), function(p) {
        temporal_subjects$person_means[p, ]
      }) # Individual node means.
    ),
    # Temporal VAR(1) coefficients (autoregressive and cross-lagged effects).
    beta = make_parameter_output(
      mean = array(population$beta_fixed, c(nnode, nnode, lag)),            # Population temporal coefficients.
      SD = array(population$intercept_temporal_re_sd_current[nnode + seq_len(n_temporal)],         # Between-person sd of each coefficient.
                 c(nnode, nnode, lag)),
      subject = lapply(seq_len(nperson), function(p) {
        array(temporal_subjects$person_betas[, , p], c(nnode, nnode, lag))                # Individual temporal coefficients.
      })
    ),
    # Covariance matrix of the person-specific node means (random intercepts).
    omega_mu = make_covariance_output(
      cov = make_parameter_output(
        mean = population$omega[seq_len(nnode), seq_len(nnode), drop = FALSE]
      )
    ),
    # Contemporaneous innovation covariance matrices.
    theta = make_covariance_output(
      cov = make_parameter_output(
        mean = population$theta_fixed,                                            # Population innovation covariance matrix.
        subject = lapply(seq_len(nperson), function(p) contemporaneous_subjects$theta[, , p])   # Individual innovation covariances.
      )
    ),
    # Joint covariance matrix of all random intercepts and temporal random slopes.
    omega = make_covariance_output(
      cov = make_parameter_output(mean = population$omega)
    ),
    # Joint correlation matrix of all random intercepts and temporal random slopes.
    omega_cor = make_covariance_output(
      cov = make_parameter_output(mean = population$omega_cor)
    )
  )

  # Store theoretical contemporaneous regression means and random-effect SDs
  # for orthogonal_re, where these are explicitly derived from the generating
  # parameters. Other methods do not use this random-effects parameterization.
  if (contemporaneous == "orthogonal_re") {

    model_parameters$gamma_theta <- make_parameter_output(
      mean = contemporaneous_subjects$gamma_theta_fixed_true,     # Population contemporaneous regression coefficients.
      SD = contemporaneous_subjects$gamma_theta_sd_true,          # Between-person sd of each contemporaneous coefficient.
      subject = contemporaneous_subjects$gamma_theta_subject_true # Individual contemporaneous regression matrices.
    )
  }



  # Calculate the population precision matrix of the contemporaneous innovations.
  kappa <- solve(population$theta_fixed)
  # Reshape the population temporal coefficients into an nnode x nnode matrix.
  beta <- matrix(population$beta_fixed, nnode, nnode)


  # Calculate the population contemporaneous partial-correlation network.
  fixed_pcc <- inv_to_net(kappa)

  # Calculate the between-person partial-correlation network from random-intercept covariances.
  omega_mu <- population$omega[seq_len(nnode), seq_len(nnode), drop = FALSE]
  if (any(diag(omega_mu) == 0)) {
    between <- omega_mu
    between[] <- NA_real_
  } else {
    between <- inv_to_net(solve(omega_mu))
  }

  # Extract the population correlation matrix of the temporal random slopes.
  random_beta_cor <- population$omega_cor[
    nnode + seq_len(n_temporal),
    nnode + seq_len(n_temporal),
    drop = FALSE
  ]


  # Store population networks in additional formats for convenient access.
  population_networks <- list(
    fixed_kappa = kappa,               # Population precision matrix of the contemporaneous innovations.
    fixed_pcc = fixed_pcc,              # Population contemporaneous partial-correlation network.
    fixed_beta = beta,                 # Population temporal network with unstandardized VAR coefficients.
    fixed_pdc = calculate_pdc(beta, population$theta_fixed, kappa), # Population temporal network expressed as partial directed correlations.
    between = between                 # Population between-person partial-correlation network.

  )

  # Store requested sparsities, realized sparsities, final random-effect SDs, and shrinking counts.
  model_parameters$generation_info <- list(
    # Number of rejected population temporal matrices before a stationary draw.
    n_beta_fixed_redraws = population$n_beta_fixed_redraws,
    # Sparsity proportions specified as simulation arguments.
    requested_sparsity = list(
      temporal = if (is.null(beta_fixed_input)) temporal_sparsity else NA_real_,
      contemporaneous = if (is.null(theta_fixed_input)) contemporaneous_sparsity else NA_real_,
      between = if (is.null(between_cov_input)) between_sparsity else NA_real_,
      random_effect_correlations = temporal_re_correlation_sparsity
    ),
    # Proportions of zero coefficients or correlations in the generated population matrices.
    realized_sparsity = list(
      temporal = 1 - network_density(beta, directed = TRUE, include_diag = TRUE,
                                     tol = sqrt(.Machine$double.eps)),
      contemporaneous = if (contemporaneous == "random_ggm") {
        mean(vapply(seq_len(nperson), function(p) {
          pcor_p <- inv_to_net(solve(matrix(contemporaneous_subjects$theta[, , p], nnode, nnode)))
          1 - network_density(pcor_p, tol = sqrt(.Machine$double.eps))
        }, numeric(1)))
      } else {1 - network_density(fixed_pcc, tol = sqrt(.Machine$double.eps))},
      between = 1 - network_density(between, tol = sqrt(.Machine$double.eps)),
      random_effect_correlations = 1 - network_density(random_beta_cor,
                                                       tol = sqrt(.Machine$double.eps))
    ),
    # Final population SDs controlling between-person variability in model parameters.
    random_effect_sd = list(
      means = population$intercept_temporal_re_sd_current[seq_len(nnode)],
      temporal = array(
        population$intercept_temporal_re_sd_current[nnode + seq_len(n_temporal)],
        c(nnode, nnode, lag)
      ),
      contemporaneous_pcor =
        population$contemporaneous_pcor_re_sd_current),
    # Count temporal and contemporaneous shrinking attempts separately.
    shrinkage = list(
      n_temporal_shrink = n_temporal_shrink,
      n_contemporaneous_shrink = n_contemporaneous_shrink)
  )


  # Assemble the simulation output without time-series data.
  simulation_output <- c(
    population_networks,
    list(
      vars = paste0("V", seq_len(nnode)),
      idvar = "ID",
      lag = lag,
      model_parameters = model_parameters
    )
  )

  # Return population and individual parameters without simulating data.
  if (only_nets) {
    return(simulation_output)
  }

  # Add simulated time-series data.
  simulation_output <- c(
    list(data = data),
    simulation_output
  )

  # Match the original helper's special output format for a single person.
  if (nperson == 1L) {
    simulation_output$data <- simulation_output$data[, -ncol(simulation_output$data), drop = FALSE]
    simulation_output$between <- NULL
    attr(simulation_output, "GVARsim") <- TRUE
    class(simulation_output) <- c("list", "GVARsim")
  } else {
    class(simulation_output) <- "mlGVARsim"
    attr(simulation_output, "mlGVARsim") <- "mlGVARsim"
  }

  # Return results
  simulation_output
}




#' Set up population parameters for a multilevel graphical VAR(1) model
#'
#' Constructs the joint random-effects covariance matrix, population node
#' means, temporal coefficients, and contemporaneous network parameters.
#' Randomly generated temporal coefficient matrices are scaled to a spectral
#' radius of at most 0.8. A supplied temporal coefficient matrix must have
#' a spectral radius below 1.
#'
#' @noRd
setup_population <- function(
    nnode, n_temporal, n_random_effects,
    random_intercept_sd_range, temporal_beta_re_sd_range,
    between_sparsity, temporal_re_correlation_sparsity,
    theta_var, contemporaneous, contemporaneous_sparsity,
    contemporaneous_pcor_re_sd_range,
    fixed_mean_sd, temporal_sparsity, temporal_prop_positive,
    beta_fixed = NULL, between_cov = NULL, theta_fixed = NULL) {

  # Random-intercept correlations define the between-person network.
  if (is.null(between_cov)) {
    omega_mu_cor <- generate_sparse_ggm_cor(nnode, between_sparsity)
  } else {
    omega_mu_cor <- stats::cov2cor(between_cov)
  }

  # Generate correlations among temporal random slopes.
  omega_beta_cor <- generate_sparse_temporal_re_cor(
    n_temporal, temporal_re_correlation_sparsity
  )

  # Random intercepts and temporal random slopes are uncorrelated.
  omega_cor <- rbind(
    cbind(omega_mu_cor, matrix(0, nnode, n_temporal)),
    cbind(matrix(0, n_temporal, nnode), omega_beta_cor)
  )

  # Draw all random-effect SDs together when both sets are generated.
  # A supplied between-person covariance determines the intercept SDs.
  if (is.null(between_cov)) {
    intercept_temporal_re_sd_initial <- stats::runif(
      n_random_effects,
      min = c(
        rep(random_intercept_sd_range[1], nnode),
        rep(temporal_beta_re_sd_range[1], n_temporal)
      ),
      max = c(
        rep(random_intercept_sd_range[2], nnode),
        rep(temporal_beta_re_sd_range[2], n_temporal)
      )
    )
  } else {
    intercept_temporal_re_sd_initial <- c(
      sqrt(diag(between_cov)),
      stats::runif(
        n_temporal,
        min = temporal_beta_re_sd_range[1],
        max = temporal_beta_re_sd_range[2]
      )
    )
  }

  # Keep the initial SDs for reference if temporal shrinkage is needed later.
  intercept_temporal_re_sd_current <- intercept_temporal_re_sd_initial
  omega <- cor_to_cov(omega_cor, intercept_temporal_re_sd_current)

  if (!is.null(between_cov)) {
    omega[seq_len(nnode), seq_len(nnode)] <- between_cov
  }

  # For random_ggm, theta_fixed is obtained later by averaging the
  # person-specific innovation covariance matrices.
  if (is.null(theta_fixed) && contemporaneous != "random_ggm") {
    theta_cor <- generate_sparse_ggm_cor(nnode, contemporaneous_sparsity)
    theta_fixed <- cor_to_cov(theta_cor, sqrt(theta_var))
  }

  # The following quantities are used only for orthogonal_re.
  k_fixed <- NULL
  k_diag <- NULL
  p_fixed <- NULL
  contemporaneous_edge_index <- NULL
  contemporaneous_edge_index_lower <- NULL
  n_contemporaneous_edges <- NULL
  contemporaneous_pcor_re_sd_current <- NULL

  if (contemporaneous == "orthogonal_re") {
    k_fixed <- solve(theta_fixed)
    k_fixed <- (k_fixed + t(k_fixed)) / 2
    k_diag <- diag(k_fixed)

    p_fixed <- inv_to_net(k_fixed)
    diag(p_fixed) <- 1

    # Include every possible edge, including edges with zero population weight.
    contemporaneous_edge_index <- which(upper.tri(p_fixed), arr.ind = TRUE)
    contemporaneous_edge_index_lower <- cbind(
      contemporaneous_edge_index[, 2],
      contemporaneous_edge_index[, 1]
    )
    n_contemporaneous_edges <- nrow(contemporaneous_edge_index)

    contemporaneous_pcor_re_sd_current <- stats::runif(
      n_contemporaneous_edges,
      min = contemporaneous_pcor_re_sd_range[1],
      max = contemporaneous_pcor_re_sd_range[2]
    )

    if (n_contemporaneous_edges > 0L) {
      names(contemporaneous_pcor_re_sd_current) <- paste0(
        "edge_",
        contemporaneous_edge_index[, 1],
        "_",
        contemporaneous_edge_index[, 2]
      )
    }
  }

  # Draw one population mean per node before generating temporal coefficients.
  mu_fixed <- stats::rnorm(nnode, mean = 0, sd = fixed_mean_sd)

  # Later simulation steps use the temporal coefficients as a vector.
  n_beta_fixed_redraws <- 0L

  if (!is.null(beta_fixed)) {
    radius <- max(Mod(eigen(beta_fixed, only.values = TRUE)$values))

    if (radius >= 1) {
      stop("Supplied beta_fixed is not stationary (spectral radius >= 1).")
    }

    beta_fixed <- as.vector(beta_fixed)
  } else {
    max_beta_fixed_attempts <- 10000L
    n_nonzero <- n_temporal - round(n_temporal * temporal_sparsity)

    for (attempt in seq_len(max_beta_fixed_attempts)) {
      beta_fixed <- numeric(n_temporal)

      if (n_nonzero > 0L) {
        nonzero_index <- sample.int(n_temporal, n_nonzero)
        signs <- sample(
          c(-1, 1), n_nonzero, replace = TRUE,
          prob = c(1 - temporal_prop_positive, temporal_prop_positive)
        )
        beta_fixed[nonzero_index] <- signs * abs(stats::rnorm(n_nonzero))
      }

      # Scale all coefficients by the same factor to preserve their relative
      # magnitudes, signs, and zero pattern.
      beta_matrix <- matrix(beta_fixed, nnode, nnode)
      radius <- max(Mod(eigen(beta_matrix, only.values = TRUE)$values))

      if (radius > 0.8) {
        beta_fixed <- beta_fixed * (0.8 / radius)
        beta_matrix <- matrix(beta_fixed, nnode, nnode)
        radius <- max(Mod(eigen(beta_matrix, only.values = TRUE)$values))
      }

      if (radius < 1) break

      n_beta_fixed_redraws <- n_beta_fixed_redraws + 1L
    }

    if (radius >= 1) {
      stop(
        "Could not generate a stationary population temporal matrix after ",
        max_beta_fixed_attempts, " attempts. Check temporal_sparsity and nnode."
      )
    }
  }

  list(
    mu_fixed = mu_fixed,
    beta_fixed = beta_fixed,
    n_beta_fixed_redraws = n_beta_fixed_redraws,
    omega_mu_cor = omega_mu_cor,
    omega_beta_cor = omega_beta_cor,
    omega_cor = omega_cor,
    intercept_temporal_re_sd_initial = intercept_temporal_re_sd_initial,
    intercept_temporal_re_sd_current = intercept_temporal_re_sd_current,
    omega = omega,
    theta_fixed = theta_fixed,
    k_fixed = k_fixed,
    k_diag = k_diag,
    p_fixed = p_fixed,
    contemporaneous_edge_index = contemporaneous_edge_index,
    contemporaneous_edge_index_lower = contemporaneous_edge_index_lower,
    n_contemporaneous_edges = n_contemporaneous_edges,
    contemporaneous_pcor_re_sd_current = contemporaneous_pcor_re_sd_current
  )
}


#' Generate person-specific innovation covariance matrices
#'
#' For orthogonal random effects, draw one set of networks for all persons.
#' If a precision matrix fails the eigenvalue check, return `valid = FALSE`
#' so that `mlvar_simulate()` can reduce the edge SDs and draw again.
#' @noRd
generate_contemporaneous_subjects <- function(
    population, nperson, nnode,
    contemporaneous, contemporaneous_sparsity,
    theta_var
) {
  # Only random_ggm generates theta_fixed here; other methods use the
  # covariance matrix already stored in the population object.
  theta_fixed <- population$theta_fixed

  # Defaults for methods that do not use independent contemporaneous edge REs.
  contemporaneous_re <- NULL
  gamma_theta_fixed_true <- NULL
  gamma_theta_sd_true <- NULL
  gamma_theta_subject_true <- NULL

  if (contemporaneous == "wishart") {
    # Set Wishart degrees of freedom locally; no other method uses them.
    df_theta <- nnode * 2

    # Generate a person-specific innovation covariance matrix for
    # each individual using Wishart distribution.
    # Dividing theta_fixed by df_theta ensures that the expected
    # value of the Wishart distribution equals theta_fixed.
    theta <- stats::rWishart(
      nperson,
      df_theta,
      theta_fixed / df_theta
    )

  } else if (contemporaneous == "random_ggm") {
    # Draw one innovation covariance matrix per person.
    innovation_sds <- sqrt(theta_var)
    theta <- array(NA_real_, dim = c(nnode, nnode, nperson))

    for (i in seq_len(nperson)) {
      theta_cor <- generate_sparse_ggm_cor(nnode, contemporaneous_sparsity)
      theta[, , i] <- cor_to_cov(theta_cor, innovation_sds)
    }

    # The population covariance is the mean of the individual covariances.
    theta_fixed <- apply(theta, c(1L, 2L), mean)

  } else if (contemporaneous == "orthogonal_re") {
    # Draw one batch of person-specific partial-correlation networks.
    # If any precision matrix fails, report the failure to the parent.

    # Draw independent standard-normal deviations for each person
    # and each unique contemporaneous edge.
    # Rows represent persons, columns represent edges.
    # Independence across columns implies uncorrelated random effects
    # between different contemporaneous edges (orthogonal_re).
    # At this stage, all deviations have mean 0 and sd 1.
    contemporaneous_re <- matrix(
      rnorm(
        nperson * population$n_contemporaneous_edges
      ),
      nrow = nperson,
      ncol = population$n_contemporaneous_edges
    )

    # Multiply the deviations for each edge by its corresponding
    # population-level random-effect sd.
    # Applies the multiplication column-wise,
    # so all persons share the same sd for a given edge.
    # The resulting values represent person-specific deviations
    # from the population partial correlations.
    contemporaneous_re <- sweep(
      contemporaneous_re,
      2,
      population$contemporaneous_pcor_re_sd_current,
      "*"
    )

    # Assign edge labels to the columns of the random-effects matrix.
    colnames(contemporaneous_re) <-
      names(population$contemporaneous_pcor_re_sd_current)

    # Initialize a list to store the person-specific innovation covariance matrices.
    theta_list <- vector(
      "list",
      nperson
    )

    # Initialize a list to store the person-specific matrices
    # of contemporaneous nodewise regression coefficients.
    gamma_theta_subject_true <- vector(
      "list",
      nperson
    )

    # Stop this batch immediately if any precision matrix is invalid.
    # The parent function then reduces the edge SDs and retries.

    for (i in seq_len(nperson)) {

      # Start from the population network, including its zero edges.
      p_i <- population$p_fixed

      # Add person i's deviation to every unique upper-triangle edge.
      # a zero population edge can therefore become nonzero for a person.
      p_i[population$contemporaneous_edge_index] <-
        p_i[population$contemporaneous_edge_index] +
        contemporaneous_re[i, ]

      # Copy the updated edges to the lower triangle to maintain
      # symmetry, and retain ones on the diagonal by convention
      p_i[population$contemporaneous_edge_index_lower] <-
        p_i[population$contemporaneous_edge_index]
      diag(p_i) <- 1

      # For j != k: K_i[j,k] = -P_i[j,k] * sqrt(K_diag[j]*K_diag[k]).
      # The outer product of the fixed precision diagonal gives pairwise products.
      # This sets the off-diagonal
      # precision entries corresponding to the individual's pcor edges.
      k_i <- -p_i * sqrt(population$k_diag %o% population$k_diag)

      # Replace the diagonal from the preceding operation with the
      # positive diagonal of the population precision matrix; it is fixed.
      diag(k_i) <- population$k_diag

      # Check positive definiteness
      min_eigenvalue <- min(
        eigen(
          k_i,
          symmetric = TRUE,
          only.values = TRUE
        )$values
      )

      # Non-finite or too-small eigenvalue: reject the ENTIRE batch
      # and let the parent function decide how to respond.
      if (
        !is.finite(min_eigenvalue) ||
        min_eigenvalue <= sqrt(.Machine$double.eps)
      ) {
        return(list(valid = FALSE))
      }

      # Invert k_i to obtain Theta_i, the covariance matrix of the VAR innovations.
      # Store it in position i of the list for later combination.

      theta_list[[i]] <- solve(
        k_i
      )

      # For a GGM, nodewise regression coefficients can be computed
      # from the precision matrix: gamma_i = I - d^(-1) k_i,
      # where d contains the diagonal entries of k_i (= population$k_diag).
      # These coefficients are generally asymmetric, unlike pcors.
      d <- diag(population$k_diag, nnode, ncol = nnode)
      # solve(d, k_i) evaluates d^(-1) %*% k_i.
      gamma_i <- diag(nnode) - solve(d, k_i)
      # Remove self-regression terms; retain only cross-node effects.
      diag(gamma_i) <- 0
      # Save the validated person's nodewise coefficient matrix.
      gamma_theta_subject_true[[i]] <- gamma_i
    }

    # After all person-specific innovation covariances pass:
    # Combine person-specific residual covariance matrices
    # Stack nnode x nnode matrices along the third dimension (base r).
    theta <- array(
      unlist(theta_list, use.names = FALSE),
      dim = c(nnode, nnode, nperson)
    )

    # Nodewise contemporaneous regression parameters for the population.
    d <- diag(population$k_diag, nnode, nnode)
    gamma_theta_fixed_true <- diag(nnode) - solve(d, population$k_fixed)
    # Again no self loops in the contemporaneous part
    diag(gamma_theta_fixed_true) <- 0

    # Initialize a matrix to store the population SDs of the
    # contemporaneous nodewise regression coefficients.
    gamma_theta_sd_true <- matrix(0, nnode, nnode)

    # Convert the random-effect SDs from partial correlations
    # to the corresponding nodewise regression coefficient scale.
    for (e in seq_len(population$n_contemporaneous_edges)) {

      # Identify the two nodes connected by the current edge.
      j <- population$contemporaneous_edge_index[e, 1]
      k <- population$contemporaneous_edge_index[e, 2]

      # Obtain the random-effect sd on the partial correlation scale.
      pcor_sd <- population$contemporaneous_pcor_re_sd_current[e]

      # Convert to the sd of the regression coefficient
      # predicting node j from node k and node k to j (not necessarily equal)
      gamma_theta_sd_true[j, k] <-
        pcor_sd * sqrt(population$k_diag[k] / population$k_diag[j])
      gamma_theta_sd_true[k, j] <-
        pcor_sd * sqrt(population$k_diag[j] / population$k_diag[k])
    }

  } else {
    # fixed: all persons share exactly the same innovation covariance matrix.
    theta <- array(
      rep(theta_fixed, nperson),
      dim = c(nnode, nnode, nperson)
    )
  }

  list(
    valid = TRUE,
    theta_fixed = theta_fixed,
    theta = theta,
    contemporaneous_re = contemporaneous_re,
    gamma_theta_fixed_true = gamma_theta_fixed_true,
    gamma_theta_sd_true = gamma_theta_sd_true,
    gamma_theta_subject_true = gamma_theta_subject_true
  )
}


#' Draw person-specific means and temporal VAR(1) coefficients
#'
#' Draws one parameter vector per person from the joint random-effects
#' distribution and checks whether all person-specific temporal matrices
#' have a spectral radius below 1.
#' @noRd
generate_temporal_subjects <- function(
    nperson, nnode, mu_fixed, beta_fixed, omega) {

  # Each row contains one person's node means followed by temporal coefficients.
  person_parameters <- mvtnorm::rmvnorm(
    n = nperson,
    mean = c(mu_fixed, beta_fixed),
    sigma = omega
  )

  person_means <- person_parameters[, seq_len(nnode), drop = FALSE]

  # Store each person's coefficients as one matrix in the third dimension.
  beta_draws <- person_parameters[, -seq_len(nnode), drop = FALSE]
  person_betas <- array(t(beta_draws), dim = c(nnode, nnode, nperson))

  # Stop checking once a nonstationary matrix is found.
  stationary <- TRUE
  for (i in seq_len(nperson)) {
    beta_i <- matrix(person_betas[, , i], nrow = nnode, ncol = nnode)
    radius <- max(Mod(eigen(beta_i, only.values = TRUE)$values))

    if (radius >= 1) {
      stationary <- FALSE
      break
    }
  }

  list(
    person_means = person_means,
    person_betas = person_betas,
    stationary = stationary
  )
}

#' Generate a correlation matrix with a sparse precision matrix
#'
#' Sparsity is the proportion of zero edges among unique off-diagonal
#' entries of the precision matrix. The number of nonzero edges is
#' rounded to the nearest integer.
#'
#' @noRd
generate_sparse_ggm_cor <- function(nnode, sparsity) {
  # Count each undirected edge once.
  precision <- matrix(0, nrow = nnode, ncol = nnode)
  upper <- upper.tri(precision)
  n_edges <- sum(upper)
  n_nonzero <- round((1 - sparsity) * n_edges)

  # Also covers the single-node case.
  if (n_nonzero == 0L) {
    return(diag(nnode))
  }

  # Select edges and draw signed weights for the upper triangle.
  weights <- numeric(n_edges)
  selected <- sample.int(n_edges, n_nonzero)
  signs <- sample(c(-1, 1), n_nonzero, replace = TRUE)
  weights[selected] <- signs * stats::runif(n_nonzero, min = 0.5, max = 1)

  precision[upper] <- weights
  precision <- precision + t(precision)

  # Strict diagonal dominance ensures positive definiteness.
  # Isolated nodes receive unit precision.
  diagonal <- 1.5 * rowSums(abs(precision))
  diagonal[diagonal == 0] <- 1
  diag(precision) <- diagonal

  # Convert to a covariance matrix, then standardize to unit variances.
  covariance <- chol2inv(chol(precision))
  stats::cov2cor(covariance)
}

#' Generate a sparse correlation matrix for temporal random slopes
#'
#' Samples marginal correlations for selected pairs of random slopes.
#' Scales the off-diagonal entries uniformly when their matrix exceeds
#' the specified spectral radius, preserving the sampled zero pattern.
#' @noRd
generate_sparse_temporal_re_cor <- function(
    n_param, sparsity, min_abs = 0.05, max_abs = 0.30,
    prop_positive = 0.50, max_radius = 0.90) {

  if (n_param == 1L || sparsity == 1) {
    return(diag(n_param))
  }

  # Count each pair of random slopes once.
  off_diag <- matrix(0, n_param, n_param)
  pairs <- which(upper.tri(off_diag), arr.ind = TRUE)
  n_pairs <- nrow(pairs)
  n_nonzero <- n_pairs - round(n_pairs * sparsity)

  if (n_nonzero == 0L) {
    return(diag(n_param))
  }

  # Draw the locations, signs, and magnitudes of the nonzero correlations.
  selected <- sample.int(n_pairs, n_nonzero)
  signs <- sample(
    c(-1, 1), n_nonzero, replace = TRUE,
    prob = c(1 - prop_positive, prop_positive)
  )
  values <- signs * stats::runif(n_nonzero, min = min_abs, max = max_abs)

  off_diag[pairs[selected, , drop = FALSE]] <- values
  off_diag <- off_diag + t(off_diag)

  # Bound the spectral radius before adding the unit diagonal.
  radius <- max(abs(eigen(off_diag, symmetric = TRUE, only.values = TRUE)$values))
  if (radius > max_radius) {
    off_diag <- off_diag * (max_radius / radius)
  }

  correlation <- diag(n_param) + off_diag
  if (min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Could not construct a positive-definite random-effect correlation matrix.",
         call. = FALSE)
  }

  correlation
}

#' Assemble population and person-specific covariance representations
#'
#' Returns covariance, correlation, precision, and partial correlation matrices.
#' For singular covariance matrices, precision and partial correlations are
#' returned as missing. Marginal correlations are computed for variables
#' with positive variance.
#' @noRd
make_covariance_output <- function(cov) {
  # The first matrix represents the population, followed by one per person.
  matrices <- c(list(cov$mean), cov$subject)

  correlations <- vector("list", length(matrices))
  precisions <- vector("list", length(matrices))
  partial_correlations <- vector("list", length(matrices))

  for (i in seq_along(matrices)) {
    x <- matrices[[i]]

    # Preserve dimensions and node names for undefined results.
    empty <- x
    empty[] <- NA_real_
    correlations[[i]] <- empty
    precisions[[i]] <- empty
    partial_correlations[[i]] <- empty

    if (anyNA(x)) next

    # Zero-variance variables have undefined marginal correlations.
    varying <- diag(x) > 0
    if (any(varying)) {
      correlations[[i]][varying, varying] <-
        stats::cov2cor(x[varying, varying, drop = FALSE])
    }

    # A zero variance makes the full covariance matrix singular.
    if (!all(varying)) next

    # Linear dependence can also occur when all variances are positive.
    precision <- try(solve(x), silent = TRUE)
    if (inherits(precision, "try-error")) next

    precisions[[i]] <- precision
    pcor <- inv_to_net(precision)
    diag(pcor) <- 1
    partial_correlations[[i]] <- pcor
  }

  # Restore the population and subject structure for each representation.
  out <- list(cov = cov)
  derived <- list(
    cor = correlations,
    prec = precisions,
    pcor = partial_correlations
  )

  for (name in names(derived)) {
    values <- derived[[name]]
    subjects <- if (is.null(cov$subject)) NULL else values[-1L]
    out[[name]] <- make_parameter_output(mean = values[[1L]], subject = subjects)
  }

  out
}


#' Collect population values, between-person SDs, and person-specific values
#' @noRd
make_parameter_output <- function(mean, SD = NULL, subject = NULL) {
  nperson <- length(subject)

  if (is.null(SD) && nperson > 1L) {
    subject_mean <- Reduce(`+`, subject) / nperson

    # Accumulate squared deviations elementwise to preserve the original shape.
    sum_sq <- (subject[[1L]] - subject_mean)^2
    for (i in 2:nperson) {
      sum_sq <- sum_sq + (subject[[i]] - subject_mean)^2
    }

    SD <- sqrt(sum_sq / (nperson - 1L))
  }

  list(mean = mean, SD = SD, subject = subject)
}


#' Simulate one person's VAR(1) time series
#'
#' Starts at `init` and returns `ntime` observations after discarding
#' the first `burnin` rows. Columns correspond to nodes.
#' @noRd
simulate_var1 <- function(beta, means, theta, ntime,
                          init = means, burnin = 100L) {
  nnode <- nrow(beta)
  n_total <- ntime + burnin

  series <- matrix(NA_real_, nrow = n_total, ncol = nnode)
  series[1L, ] <- init

  if (n_total > 1L) {
    innovations <- mvtnorm::rmvnorm(
      n = n_total - 1L, mean = numeric(nnode), sigma = theta
    )

    for (i in seq.int(2L, n_total)) {
      # Apply temporal effects to deviations from the person's long-run means.
      previous <- series[i - 1L, ] - means
      series[i, ] <- means + as.vector(beta %*% previous) + innovations[i - 1L, ]
    }
  }

  series[burnin + seq_len(ntime), , drop = FALSE]
}
