## util
compute_icate <- function(fit, test = NULL){
  if(class(fit)[1] == 'bart'){
    post <- dbarts::extract(fit, 'ev', 'test')
  }
  
  if(class(fit)[1] == 'bcfmodel'){
    return(t(fit$tau_hat_train)) # transpose becuase stochtree index is different

  }
  
  if(class(fit)[1] == 'stanreg'){
    post <- rstanarm::posterior_epred(fit, newdata = test)
  }
  
  if(class(fit)[1] == 'stan4bartFit'){
    post <- t(stan4bart:::extract.stan4bartFit(fit, 'ev', 'test'))
  }
  
   post[, 1:(nrow(test)/2)] - post[, ((nrow(test)/2)+1):nrow(test)]
  
  
  
}

compute_subgroup_effects <- function(n_groups = 8, cate_truth = cate_truth, icate = icate, dat = dat, model){
  absolute_cates <- lapply(1:n_groups, function(j){
    cate <- apply(icate[, dat$x == j], 1, mean)
    data.frame(truth = cate_truth$coef[[j]], 
               est = mean(cate), 
               lci = quantile(cate, .025), 
               uci = quantile(cate, .975), 
               sd_y = sd(dat$y),
               model,
               type = 'absolute',
               row.names = NULL)
  })
  
  absolute_cates <- do.call('rbind', absolute_cates)
  # compute all possible contrast effects
  combinations <- as.data.frame(t(combn(n_groups, 2)))
  # compute relative subgroup effects
  relative_cates <- lapply(1:nrow(combinations), function(j){
    g1 <- apply(icate[, dat$x == combinations$V1[j]], 1, mean)
    g2 <- apply(icate[, dat$x == combinations$V2[j]], 1, mean)
    cate <-  g1 - g2 
    data.frame(truth = cate_truth$coef[[combinations$V1[j]]] - cate_truth$coef[[combinations$V2[j]]], 
               est = mean(cate), 
               lci = quantile(cate, .025), 
               uci = quantile(cate, .975), 
               sd_y = sd(dat$y),
               model,
               type = 'relative',
               row.names = NULL)
  })
  
  relative_cates <- do.call('rbind', relative_cates)
  rbind(absolute_cates, relative_cates)
}
