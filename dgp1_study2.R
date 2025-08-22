source('utils.R')
iter <- Sys.getenv("SLURM_ARRAY_TASK_ID")

# participants in the original IHDP study
n <- 985 

# site counts from IHDP study for sihgts 1-8
site_counts <- c(128, 138, 138, 100, 101, 137, 131, 112)

# x is a vector of site membership values 1-8
x <- lapply(seq_along(site_counts), function(i){rep(i, site_counts[i])})
x <- unlist(x)

# create x_mat a matrix of indicator variables
x_mat <- matrix(NA, nrow = length(x), ncol = 8)
for (i in 1:8) {
  x_mat[, i] <- x == i
}

x_mat <- apply(x_mat, 2, as.double)
# grid of subgroup effect sizes
d <- seq(0, 1.5, .1) 

run <- lapply(1:length(d), function(i){
  set.seed(2)
  # mean structure for potenital outcomes
  y0hat <- rep(0, n)
  
  # note everything is 0 expet for group 8 which is set to d
  beta <- rnorm(8, 0, d[i])
  # assign subgroup effect for group 8
  y1hat <- as.vector(x_mat%*%beta)
  
  # compute true subgroup effects 
  cate_truth <- y1hat - y0hat
  cate_truth <- lm(cate_truth ~ 0 + as.factor(x))
  
  # draw random z and error values, create Y(1) and Y(0)
  set.seed(iter)
  z <- rbinom(n, 1, .5)
  y0 <- y0hat + rnorm(n)
  y1 <- y1hat + rnorm(n)
  
  # y is the observed potential outcome
  y <- ifelse(z == 1, y1, y0)
  
  # put it all in a data.frame for book keeping
  dat <- data.frame(y, z, x)
  
  # create a test matrix, we'll use this to compute CATE estimates for stan/BART
  test_df <- rbind(dat, dat)
  test_df$z <- c(rep(1, nrow(dat)), rep(0, nrow(dat)))
  
  # run models
  fit1 <- rstanarm::stan_glmer(y ~ z + (1 + z|x), data = dat)
  fit2 <- dbarts::bart2(y ~ z + as.factor(x), data = dat, test = test_df)
  fit3 <- stochtree::bcf(X_train = x_mat, Z_train = z, y_train = y, num_mcmc = 1000)
  fit4 <- rstanarm::stan_glm(y ~ z*as.factor(x), data = dat)
  
  # extract icate (we'll need this for SVD)
  icate_fit1 <- compute_icate(fit1, test = test_df)
  icate_fit2 <- compute_icate(fit2, test = test_df)
  icate_fit3 <- compute_icate(fit3)
  icate_fit4 <- compute_icate(fit4, test = test_df)
  
  # save as a list for later
  icate <- list(mlm_icate = icate_fit1, 
                BART_icate = icate_fit2, 
                bcf_icate = icate_fit3, 
                lm_icate = icate_fit4)
  
  # extract statistics
  subgroup_effects <- rbind(
    compute_subgroup_effects(
      n_groups = 8,
      cate_truth = cate_truth,
      icate = icate_fit1,
      dat = dat,
      model = 'mlm'
    ),
    
    compute_subgroup_effects(
      n_groups = 8,
      cate_truth = cate_truth,
      icate = icate_fit2,
      dat = dat,
      model = 'BART'
    ),
    
    compute_subgroup_effects(
      n_groups = 8,
      cate_truth = cate_truth,
      icate = icate_fit3,
      dat = dat,
      model = 'bcf'
    ),
    
    compute_subgroup_effects(
      n_groups = 8,
      cate_truth = cate_truth,
      icate = icate_fit4,
      dat = dat,
      model = 'lm'
    )
  )
  
  # keep track of iteration and d for book keeping 
  subgroup_effects$iter <- iter
  subgroup_effects$d <- d[i]
  
  list(subgroup = subgroup_effects, icate = icate)
  
})


export_subgroups <- do.call('rbind', lapply(run, '[[', 1))
export_icates <- lapply(run, '[[', 2)

path <- paste0('./results/subgroups_dgp1_study2/results', iter, '.csv')
readr::write_csv(export_subgroups, path)

path <- paste0('./results/icates_dgp1_study2/results', iter, '.rds')
readr::write_rds(export_icates, path)
