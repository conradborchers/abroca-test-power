# Step 1
generate_data <- function(n=1000, auc=0.8, ratio_pos_case=0.5) {
  # Numbers based on proof in
  # SALGADO, Jesús F.. Transforming the Area under the Normal Curve (AUC) 
  #into Cohen’s d, Pearson’s r pb , Odds-Ratio, and Natural Log Odds-Ratio: 
  #Two Conversion Tables. The European Journal of Psychology Applied to 
  #Legal Context [online]. 2018, vol.10, n.1, pp.35-47. ISSN 1989-4007. 
  #http://dx.doi.org/10.5093/ejpalc2018a5
  # See also: https://stats.stackexchange.com/questions/422926/generate-synthetic-data-given-auc
  t <- sqrt(log(1/(1-auc)**2))
  z <- t-((2.515517 + 0.802853*t + 0.0103328*t**2) / 
            (1 + 1.432788*t + 0.189269*t**2 + 0.001308*t**3))
  d <- z*sqrt(2)
  x <- c(rnorm(n*(1-ratio_pos_case), mean = 0), rnorm(n*ratio_pos_case, mean = d))
  y <- c(rep(0, n*(1-ratio_pos_case)), rep(1, n*ratio_pos_case))
  return(data.frame(x=x, y=y))
}

# Step 2
generate_predictions <- function(n_majority, n_minority, auc_majority, auc_minority, ratio_pos_case, test_set_size=0.2) {
  data0 <- generate_data(n_minority, auc_minority, ratio_pos_case)
  data0$group <- 0
  data1 <- generate_data(n_majority, auc_majority, ratio_pos_case)
  data1$group <- 1
  data <- rbind(data0, data1)
  
  train_idx <- sample(1:nrow(data), (1-test_set_size) * nrow(data))
  train_data <- data[train_idx, ]
  test_data <- data[-train_idx, ]
  
  model <- glm(y ~ x, data = train_data, family = binomial)
  pred <- predict(model, test_data, type = "response")
  
  test_data['pred'] <- pred
  return(test_data)
}

calculate_abroca <- function(test_data, randomize_groups=FALSE, return_both_aucs=FALSE) {
  if (randomize_groups)
    test_data$group <- sample(test_data$group, nrow(test_data), replace=FALSE)
  roc_all <- pROC::roc(test_data$y, test_data$pred)
  roc_0 <- pROC::roc(test_data$y[test_data$group == 0], test_data$pred[test_data$group == 0])
  roc_1 <- pROC::roc(test_data$y[test_data$group == 1], test_data$pred[test_data$group == 1])
  tpr_0_interp <- stats::approx(roc_0$specificities, roc_0$sensitivities, xout = roc_all$specificities)$y
  tpr_1_interp <- stats::approx(roc_1$specificities, roc_1$sensitivities, xout = roc_all$specificities)$y
  abroca <- pracma::trapz(roc_all$specificities, abs(tpr_0_interp - tpr_1_interp))
  if (return_both_aucs) {
    return(list(abroca=abroca, auc0=as.numeric(roc_0$auc), auc1=as.numeric(roc_1$auc)))
  }
  return(abroca)
}

simulate_abroca_lak <- function(n_majority, n_minority, auc_majority, auc_minority, ratio_pos_case, test_set_size=0.2, return_both_aucs=FALSE) {
  test_data <- generate_predictions(n_majority, n_minority, auc_majority, auc_minority, ratio_pos_case, test_set_size=0.2)
  abroca_obs <- calculate_abroca(test_data, return_both_aucs=return_both_aucs)
  return(abroca_obs)
}

# Steps 3 and 4
simulate_abroca <- function(n_majority, n_minority, auc_majority, auc_minority, ratio_pos_case, test_set_size=0.2, n_replications=1000, return_both_aucs=FALSE) {
  test_data <- generate_predictions(n_majority, n_minority, auc_majority, auc_minority, ratio_pos_case, test_set_size=0.2)
  abroca_obs <- calculate_abroca(test_data, return_both_auc=return_both_aucs)
  abroca_null <- replicate(n_replications, {
    tryCatch({
      suppressMessages({calculate_abroca(test_data, randomize_groups=TRUE, return_both_aucs=return_both_aucs)})
    }, error=function(e){NA}) # Discuss
  })
  # Naturally two-sided as abroca cannot be positive and is two-sided per definition
  p <- mean(abroca_null>abroca_obs, na.rm=TRUE)
  return(p)
}

main <- function(n_majority, n_minority, auc_majority, auc_minority, ratio_pos_case, test_set_size=0.2, n_replications=1000) {
  ps <- replicate(n_replications, { 
    tryCatch({
      ans <- suppressMessages({simulate_abroca(n_majority, n_minority, auc_majority, auc_minority, ratio_pos_case=ratio_pos_case, test_set_size=test_set_size, n_replications=100)})
      ans 
    }, error=function(e){NA}) # Can impute 0 or remove as per the user's preference if data set is too small for having both classes per group
  })
  return(ps)
}

run_and_store_simulation <- function(ratio_minority, auc_majority, auc_minority, ratio_pos_case, test_set_size=0.2, n_replications=1000) {
  results <- list()
  for (n in seq(500, 10000, 500)) { # Defaults
    cat('Simulating for n =', n, '\n')
    ps <- main(n_majority=n*(1-ratio_minority), n_minority=n*ratio_minority, 
               auc_majority=auc_majority, auc_minority=auc_minority, 
               ratio_pos_case=ratio_pos_case, 
               test_set_size=test_set_size, n_replications=n_replications) 
    results[[as.character(n)]] <- ps 
  }
  f <- paste('ratio_minority_',ratio_minority, '_auc_majority_', auc_majority,
             '_auc_minority_', auc_minority, '_ratio_pos_case_', ratio_pos_case, '-pvals-1000-v1.rds', sep='')
  saveRDS(results, file=f)
  return(results)
}

# Evaluation

doEval <- function(f_rds='ratio_minority_0.5_auc_majority_0.8_auc_minority_0.65_ratio_pos_case_0.1-v1.rds') {
  sim <- readRDS(f_rds)
  ref <- readRDS(f_rds %>% str_replace('auc_minority_0.[0-9]{1,}', 'auc_minority_0.8'))
  threshold <- ref[max(names(ref))] %>% unlist() %>% median(na.rm=TRUE)
  join_this <- ref %>% 
    sapply(function(x){median(x, na.rm=TRUE)}) %>% 
    as.data.frame() %>% 
    rownames_to_column('n') %>%
    `colnames<-`(c('n', 'baseline'))
  p_power <- map2(sim, ref, function(x, y){return(mean(x>median(y, na.rm=TRUE), na.rm=TRUE))}) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('n') %>%
    mutate(n = n %>% str_remove('X')) %>%
    mutate(power = paste(round(V1*100, 2), '%', sep='')) %>%
    select(-V1)
  d_plot <- sim %>% 
    sapply(function(x){quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE)}) %>% 
    as.data.frame() %>% t() %>%
    as.data.frame() %>%
    rownames_to_column('n') %>%
    left_join(join_this, by='n') %>%
    left_join(p_power, by='n')
  d_plot
  d_plot %>%
    mutate(n = as.numeric(n)) %>% 
    janitor::clean_names() %>%
    ggplot(aes(n, x50_percent)) + 
    geom_point() + 
    geom_line(aes(y = baseline), color = "red") + 
    geom_text(aes(label = power, y = x97_5_percent + 0.05), size = 3) + # Adjust position as needed
    geom_errorbar(aes(ymin = x2_5_percent, ymax = x97_5_percent)) + 
    labs(
      x='Total Sample Size',
      y='ABROCA (With Red Equal AUC Baseline)',
      title=f_rds
    ) +
    theme_bw()
  
}

