#This code estimates restricted mean survival time and expected offspring under recent evidence of infection and no recent evidence of #infection. the script resamples complete male tenure histories, refits the primary cox model across 1,000 male-level bootstrap replicates, #repeats the standardized predictions, and calculates 95% percentile bootstrap confidence intervals. this analysis was run using parallel #processing on ASU’s Sol high-performance computing cluster.

# load packages
library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(future)
library(furrr)

# read the final dataset
tenure = readRDS("tenure_bootstrap_input.rds")

# confirm factor levels
tenure = tenure %>% mutate(
  infection_evidence = factor(infection_evidence, levels = c("not_detected", "detected")),
  group_size = factor(group_size, levels = c("small", "medium", "large"))
)

# fit the primary model
cm1 = coxph(Surv(start, stop, event) ~ infection_evidence + age_at_month + total_follower_males + strata(group_size), data = tenure, ties = "efron", id = tenure.id, cluster = code, x = TRUE)

# calculate rmst from a step survival curve
rmst_from_curve = function(time, surv, tau) {
  curve = tibble(time = as.numeric(time), surv = as.numeric(surv)) %>% filter(is.finite(time), is.finite(surv), time > 0, time < tau) %>% arrange(time) %>% distinct(time, .keep_all = TRUE)
  if (nrow(curve) == 0) return(tau)
  sum(diff(c(0, curve$time, tau)) * c(1, curve$surv))
}

# create empirical tenure-start profiles, giving each male equal total weight
make_prediction_profiles = function(dat) {
  dat %>%
    filter(start == 0) %>%
    distinct(code, tenure.id, .keep_all = TRUE) %>%
    transmute(code, tenure.id, starting_age = age_at_month, group_size, female_number = N.Females) %>%
    add_count(code, name = "n_profiles") %>%
    mutate(profile_weight = 1 / (n_distinct(code) * n_profiles), profile_id = sprintf("profile_%05d", row_number()))
}

# calculate standardized rmst and expected offspring under both infection conditions
compute_standardized_rmst = function(fit, dat, taus, interbirth_interval_years) {
  profiles = make_prediction_profiles(dat)
  if (nrow(profiles) == 0) stop("the bootstrap sample contains no prediction profiles")
  
  reference_age = weighted.mean(profiles$starting_age, profiles$profile_weight)
  
  reference_paths = profiles %>%
    distinct(group_size) %>%
    arrange(group_size) %>%
    mutate(reference_id = paste0("reference_", as.character(group_size)))
  
  reference_newdata = expand_grid(reference_paths, start = 0:(ceiling(max(taus)) - 1)) %>%
    mutate(
      stop = start + 1,
      event = 0,
      infection_evidence = factor("not_detected", levels = levels(dat$infection_evidence)),
      age_at_month = reference_age + start / 12,
      total_follower_males = 0,
      group_size = factor(group_size, levels = levels(dat$group_size))
    )
  
  predicted_survival = survfit(fit, newdata = reference_newdata, id = reference_id, se.fit = FALSE)
  curve_counts = as.integer(predicted_survival$strata)
  curve_order = reference_paths$reference_id
  
  if (length(curve_counts) != nrow(reference_paths) || sum(curve_counts) != length(predicted_survival$surv)) {
    stop("the predicted survival object does not contain the expected reference curves")
  }
  
  curve_names = names(predicted_survival$strata)
  
  if (!is.null(curve_names)) {
    named_order = map_chr(curve_names, function(curve_name) {
      matched_paths = curve_order[vapply(curve_order, function(path) grepl(path, curve_name, fixed = TRUE), logical(1))]
      if (length(matched_paths) == 1) matched_paths else NA_character_
    })
    if (!any(is.na(named_order))) curve_order = named_order
  }
  
  reference_curves = tibble(
    time = predicted_survival$time,
    reference_surv = as.numeric(predicted_survival$surv),
    reference_id = rep(curve_order, times = curve_counts)
  ) %>%
    left_join(reference_paths, by = "reference_id")
  
  profile_states = bind_rows(
    profiles %>% mutate(infection_state = "not_detected"),
    profiles %>% mutate(infection_state = "detected")
  ) %>%
    mutate(infection_state = factor(infection_state, levels = c("not_detected", "detected")))
  
  beta_age = unname(coef(fit)[["age_at_month"]])
  beta_infection = unname(coef(fit)[["infection_evidencedetected"]])
  
  if (!is.finite(beta_age) || !is.finite(beta_infection)) stop("the fitted model contains non-finite coefficients")
  
  reference_curves %>%
    select(time, reference_surv, group_size) %>%
    inner_join(profile_states, by = "group_size", relationship = "many-to-many") %>%
    mutate(relative_hazard = exp(beta_age * (starting_age - reference_age) + beta_infection * (infection_state == "detected")), surv = reference_surv^relative_hazard) %>%
    reframe(tau = taus, rmst = map_dbl(taus, function(tau_value) rmst_from_curve(time, surv, tau_value)), .by = c(profile_id, infection_state)) %>%
    left_join(profile_states %>% select(profile_id, infection_state, profile_weight, female_number), by = c("profile_id", "infection_state")) %>%
    mutate(offspring = rmst * female_number / (interbirth_interval_years * 12)) %>%
    summarise(rmst = weighted.mean(rmst, profile_weight), offspring = weighted.mean(offspring, profile_weight), .by = c(tau, infection_state)) %>%
    pivot_wider(names_from = infection_state, values_from = c(rmst, offspring)) %>%
    transmute(
      tau,
      rmst_not_detected,
      rmst_detected,
      lost_months = rmst_not_detected - rmst_detected,
      offspring_not_detected,
      offspring_detected,
      offspring_lost = offspring_not_detected - offspring_detected
    ) %>%
    arrange(tau)
}

# define time horizons
average_tenure_horizon = tenure %>% summarise(tenure_length_months = max(stop) - min(start), .by = tenure.id) %>% summarise(value = mean(tenure_length_months)) %>% pull(value)
plot_taus = c(12, 24, 36, 48, 60)
taus = sort(unique(c(plot_taus, average_tenure_horizon)))
interbirth_interval_years = 2.57

# calculate point estimates from the original sample
rmst_point = compute_standardized_rmst(cm1, tenure, taus, interbirth_interval_years)

# configure parallel processing using the cores assigned by slurm
workers = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
options(future.globals.maxSize = 2 * 1024^3)
plan(multisession, workers = workers)

# bootstrap complete male histories and repeat the full estimation procedure
B = 1000
male_codes = unique(tenure$code)

message("starting ", B, " bootstrap samples using ", workers, " workers")

boot_res = future_map_dfr(
  seq_len(B),
  function(b) {
    bootstrap_index = tibble(original_code = sample(male_codes, length(male_codes), replace = TRUE), bootstrap_copy = seq_along(original_code))
    
    bootstrap_data = bootstrap_index %>%
      inner_join(tenure %>% rename(original_code = code, original_tenure = tenure.id), by = "original_code", relationship = "many-to-many") %>%
      mutate(
        code = paste0(original_code, "_", bootstrap_copy),
        tenure.id = paste0(original_tenure, "_", bootstrap_copy),
        infection_evidence = factor(infection_evidence, levels = levels(tenure$infection_evidence)),
        group_size = factor(group_size, levels = levels(tenure$group_size))
      )
    
    tryCatch({
      bootstrap_fit = suppressWarnings(
        coxph(
          Surv(start, stop, event) ~ infection_evidence + age_at_month + total_follower_males + strata(group_size),
          data = bootstrap_data,
          ties = "efron",
          id = tenure.id,
          cluster = code,
          x = TRUE
        )
      )
      
      if (any(!is.finite(coef(bootstrap_fit)))) stop("non-finite model coefficient")
      
      compute_standardized_rmst(bootstrap_fit, bootstrap_data, taus, interbirth_interval_years) %>% mutate(rep = b)
    }, error = function(e) tibble())
  },
  .options = furrr_options(seed = 123, packages = c("survival", "dplyr", "tidyr", "purrr"))
)

plan(sequential)

# stop if every bootstrap sample fails
if (nrow(boot_res) == 0 || !"rep" %in% names(boot_res)) stop("all bootstrap samples failed")

# inspect bootstrap completion
bootstrap_completion = tibble(
  requested_bootstraps = B,
  successful_bootstraps = n_distinct(boot_res$rep),
  percent_successful = round(100 * n_distinct(boot_res$rep) / B, 1)
)

print(bootstrap_completion)

if (bootstrap_completion$percent_successful < 95) warning("fewer than 95% of bootstrap samples produced valid estimates")

# calculate percentile bootstrap confidence intervals
rmst_ci = boot_res %>%
  summarise(
    rmst_not_detected_lwr = quantile(rmst_not_detected, 0.025, na.rm = TRUE),
    rmst_not_detected_upr = quantile(rmst_not_detected, 0.975, na.rm = TRUE),
    rmst_detected_lwr = quantile(rmst_detected, 0.025, na.rm = TRUE),
    rmst_detected_upr = quantile(rmst_detected, 0.975, na.rm = TRUE),
    lost_months_lwr = quantile(lost_months, 0.025, na.rm = TRUE),
    lost_months_upr = quantile(lost_months, 0.975, na.rm = TRUE),
    offspring_not_detected_lwr = quantile(offspring_not_detected, 0.025, na.rm = TRUE),
    offspring_not_detected_upr = quantile(offspring_not_detected, 0.975, na.rm = TRUE),
    offspring_detected_lwr = quantile(offspring_detected, 0.025, na.rm = TRUE),
    offspring_detected_upr = quantile(offspring_detected, 0.975, na.rm = TRUE),
    offspring_lost_lwr = quantile(offspring_lost, 0.025, na.rm = TRUE),
    offspring_lost_upr = quantile(offspring_lost, 0.975, na.rm = TRUE),
    .by = tau
  )

# combine point estimates and confidence intervals
rmst_full = rmst_point %>%
  left_join(rmst_ci, by = "tau") %>%
  mutate(horizon = if_else(near(tau, average_tenure_horizon), "mean observed tenure", "additional horizon"))

# create the formatted results table
rmst_results = rmst_full %>%
  transmute(
    time_horizon = round(tau, 2),
    horizon,
    `rmst no recent evidence` = sprintf("%.2f (%.2f–%.2f)", rmst_not_detected, rmst_not_detected_lwr, rmst_not_detected_upr),
    `rmst recent evidence` = sprintf("%.2f (%.2f–%.2f)", rmst_detected, rmst_detected_lwr, rmst_detected_upr),
    `months lost` = sprintf("%.2f (%.2f–%.2f)", lost_months, lost_months_lwr, lost_months_upr),
    `offspring lost` = sprintf("%.2f (%.2f–%.2f)", offspring_lost, offspring_lost_lwr, offspring_lost_upr)
  )

print(rmst_results)

# save all results
write.csv(bootstrap_completion, "bootstrap_completion.csv", row.names = FALSE)
write.csv(rmst_results, "rmst_results.csv", row.names = FALSE)
saveRDS(boot_res, "bootstrap_replicates.rds")
saveRDS(rmst_full, "rmst_full_results.rds")
capture.output(sessionInfo(), file = "session_info.txt")

message("bootstrap analysis completed")
