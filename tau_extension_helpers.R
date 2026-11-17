####### defining shared helpers for tau extensions #######

read_numeric_csv <- function(path) {
  scan(path, what = numeric(), quiet = TRUE, sep = ",")
}

posterior_median <- function(directory, parameter, tail_n = NULL) {
  path <- file.path(directory, paste0(parameter, "_mc.csv"))
  x <- read_numeric_csv(path)
  if (!is.null(tail_n) && length(x) > tail_n) x <- tail(x, tail_n)
  unname(median(x, na.rm = TRUE))
}

haversine_matrix <- function(recipient_coords, donor_coords) {
  to_rad <- pi / 180
  lon1 <- recipient_coords[, 1] * to_rad
  lat1 <- recipient_coords[, 2] * to_rad
  lon2 <- donor_coords[, 1] * to_rad
  lat2 <- donor_coords[, 2] * to_rad
  dlon <- outer(lon1, lon2, "-")
  dlat <- outer(lat1, lat2, "-")
  a <- sin(dlat / 2)^2 + outer(cos(lat1), cos(lat2)) * sin(dlon / 2)^2
  a[a < 0] <- 0
  a[a > 1] <- 1
  6371.0088 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

read_vaccine_era <- function(repository_root) {
  directory <- file.path(repository_root, "data", "formatted", "vac")
  cases <- read.csv(
    file.path(directory, "inferred_cases_urban.csv"), check.names = FALSE
  )
  population <- read.csv(
    file.path(directory, "inferred_pop_urban.csv"), check.names = FALSE
  )
  susceptibles <- read.csv(
    file.path(directory, "susceptibles_urban.csv"), check.names = FALSE
  )
  coordinates <- read.csv(
    file.path(directory, "coordinates_urban.csv"), check.names = FALSE
  )

  drop <- intersect(c("Yeovil.RD", "City.of.London"), names(cases))
  keep <- !names(cases) %in% drop
  cases <- cases[, keep, drop = FALSE]
  population <- population[, keep, drop = FALSE]
  susceptibles <- susceptibles[, keep, drop = FALSE]
  coordinates <- coordinates[, keep, drop = FALSE]

  list(
    time = cases[[1]],
    cases = as.matrix(cases[, -1, drop = FALSE]),
    population = as.matrix(population[, -1, drop = FALSE]),
    susceptibles = as.matrix(susceptibles[, -1, drop = FALSE]),
    coordinates = t(as.matrix(coordinates[, -1, drop = FALSE])),
    names = names(cases)[-1]
  )
}

make_positive_definite <- function(x, minimum = 1e-10) {
  eig <- eigen(x, symmetric = TRUE)
  eig$values[eig$values < minimum] <- minimum
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

safe_covariance <- function(hessian, fallback_scale = 0.05) {
  covariance <- try(solve(hessian), silent = TRUE)
  if (inherits(covariance, "try-error") ||
      any(!is.finite(covariance))) {
    covariance <- diag(fallback_scale^2, nrow(hessian))
  }
  make_positive_definite(covariance)
}

laplace_draws <- function(mode, hessian, n_draw = 2500L, seed = 20201186L) {
  covariance <- safe_covariance(hessian)
  set.seed(seed)
  z <- matrix(rnorm(as.integer(n_draw) * length(mode)), ncol = length(mode))
  sweep(z %*% chol(covariance), 2, mode, "+")
}

standardize_vector <- function(x) {
  center <- mean(x, na.rm = TRUE)
  scale <- sd(x, na.rm = TRUE)
  stopifnot(is.finite(center), is.finite(scale), scale > 0)
  list(values = (x - center) / scale, center = center, scale = scale)
}

period_from_time <- function(time) {
  ifelse(
    time < 65, 1L,
    ifelse(time < 80, 2L, ifelse(time < 85, 3L, ifelse(time < 90, 4L, 5L)))
  )
}

period_labels <- c("1944-64", "1965-79", "1980-84", "1985-89", "1990-94")

score_from_lambda <- function(lambda, event) {
  probability <- -expm1(-pmax(lambda, 0))
  probability <- pmin(pmax(probability, 1e-12), 1 - 1e-12)
  data.frame(
    log_score = sum(ifelse(event, log(probability), log1p(-probability))),
    brier_score = mean((as.numeric(event) - probability)^2)
  )
}

calibration_table <- function(probability, event, model, n_bin = 10L) {
  breaks <- unique(quantile(
    probability,
    probs = seq(0, 1, length.out = n_bin + 1L),
    na.rm = TRUE
  ))
  if (length(breaks) < 3L) return(data.frame())
  group <- cut(probability, breaks = breaks, include.lowest = TRUE)
  out <- aggregate(
    cbind(predicted = probability, observed = as.numeric(event)),
    by = list(bin = group),
    FUN = mean
  )
  out$n <- as.integer(table(group)[as.character(out$bin)])
  out$model <- model
  out
}

make_spatial_regions <- function(coords, n_region = 8L, seed = 20201186L) {
  n_region <- as.integer(n_region)
  stopifnot(n_region >= 2L, n_region < nrow(coords))
  xy <- scale(coords)
  set.seed(seed)
  fit <- kmeans(xy, centers = n_region, nstart = 50, iter.max = 100)
  centers <- aggregate(
    coords,
    by = list(old_region = fit$cluster),
    FUN = mean
  )
  order_centers <- order(-centers[, 3], centers[, 2])
  remap <- integer(n_region)
  remap[centers$old_region[order_centers]] <- seq_len(n_region)
  list(
    region = remap[fit$cluster],
    centers = centers[order_centers, -1, drop = FALSE]
  )
}

prepare_tau_extension_data <- function(repository_root, n_region = 8L,
                                       seed = 20201186L) {
  vac <- read_vaccine_era(repository_root)
  vac_mcmc <- file.path(
    repository_root, "output", "mcmc inference", "vac"
  )
  recipients <- as.character(read.csv(
    file.path(vac_mcmc, "small_cities.csv"), check.names = FALSE
  )[[1]])
  recipient_index <- match(recipients, vac$names)
  stopifnot(!anyNA(recipient_index))

  fixed_parameters <- list(
    beta = posterior_median(vac_mcmc, "beta_coupling"),
    alpha_1 = posterior_median(vac_mcmc, "alpha_1"),
    alpha_2 = posterior_median(vac_mcmc, "alpha_2"),
    rho = posterior_median(vac_mcmc, "rho"),
    eta = posterior_median(vac_mcmc, "eta")
  )

  baseline_tau <- c(
    posterior_median(vac_mcmc, "tau"),
    posterior_median(vac_mcmc, "tau_2"),
    posterior_median(vac_mcmc, "tau_3"),
    posterior_median(vac_mcmc, "tau_4"),
    posterior_median(vac_mcmc, "tau_5")
  )

  recipient_coords <- vac$coordinates[recipient_index, , drop = FALSE]
  distance <- haversine_matrix(recipient_coords, vac$coordinates)
  distance[distance < 1e-8] <- Inf
  spatial_kernel <- distance^(-fixed_parameters$rho)
  donor_signal <- t(vac$cases^fixed_parameters$alpha_1)
  gravity_component <- fixed_parameters$beta *
    t(vac$population[, recipient_index, drop = FALSE]^fixed_parameters$alpha_2) *
    (spatial_kernel %*% donor_signal)
  susceptibility_factor <- t(
    pmax(
      vac$susceptibles[, recipient_index, drop = FALSE],
      1e-12
    )^fixed_parameters$eta
  )

  recipient_cases <- t(vac$cases[, recipient_index, drop = FALSE])
  n_recipient <- nrow(recipient_cases)
  n_time <- ncol(recipient_cases)
  at_risk <- recipient_cases[, seq_len(n_time - 1L), drop = FALSE] == 0
  event <- at_risk & recipient_cases[, 2:n_time, drop = FALSE] > 0

  period_time <- period_from_time(vac$time[seq_len(n_time - 1L)])
  city_all <- rep(seq_len(n_recipient), times = n_time - 1L)
  period_all <- rep(period_time, each = n_recipient)
  calendar_all <- rep(
    1900 + vac$time[seq_len(n_time - 1L)],
    each = n_recipient
  )
  susceptibility_all <- as.vector(
    susceptibility_factor[, seq_len(n_time - 1L), drop = FALSE]
  )
  gravity_all <- as.vector(
    gravity_component[, seq_len(n_time - 1L), drop = FALSE]
  )
  at_risk_all <- as.vector(at_risk)
  event_all <- as.vector(event)

  included <- at_risk_all & is.finite(susceptibility_all) &
    is.finite(gravity_all)

  london <- matrix(c(-0.1276, 51.5072), nrow = 1)
  distance_london <- as.numeric(haversine_matrix(recipient_coords, london))
  distance_standardization <- standardize_vector(distance_london)
  spatial_regions <- make_spatial_regions(
    recipient_coords, n_region = n_region, seed = seed
  )

  list(
    vac = vac,
    recipients = recipients,
    recipient_index = recipient_index,
    recipient_coords = recipient_coords,
    fixed_parameters = fixed_parameters,
    baseline_tau = baseline_tau,
    gravity_component = gravity_component,
    susceptibility_factor = susceptibility_factor,
    recipient_cases = recipient_cases,
    n_recipient = n_recipient,
    n_time = n_time,
    at_risk = at_risk,
    event = event,
    city = city_all[included],
    period = period_all[included],
    calendar = calendar_all[included],
    susceptibility = susceptibility_all[included],
    gravity = gravity_all[included],
    event_risk = event_all[included],
    distance_london = distance_london,
    distance_london_z = distance_standardization$values,
    distance_london_center = distance_standardization$center,
    distance_london_scale = distance_standardization$scale,
    region = spatial_regions$region,
    region_centers = spatial_regions$centers
  )
}
