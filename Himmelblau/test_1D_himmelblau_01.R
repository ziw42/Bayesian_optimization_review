### 1D version of the himmelblau which only studies the behavior with
### respect to the horizontal coordinate

library(tidyverse)

### define the function
himmelblau <- function(x, y)
{
  (x^2 + y - 11)^2 + (x + y^2 - 7)^2
}

### visualize the surface
df_viz <- expand.grid(x = seq(-5, 5, length.out = 301),
                      y = seq(-5, 5, length.out = 301),
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE) %>% 
  as.data.frame() %>% tibble::as_tibble()

df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = f)) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

### visualize the output in the log-scale
df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = log10(f + 1))) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

### identify the 4 modes
df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = log10_f)) +
  geom_point(data = df_viz %>% 
               mutate(f = himmelblau(x, y)) %>% 
               mutate(log_f = log10(f + 1)) %>% 
               arrange(log_f) %>% 
               slice(1:4),
             shape = 0, color = 'white', size = 3) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

### create a 1D slice at a fixed vertical position
expand.grid(x = seq(-5, 5, length.out = 301),
            y = -2.5,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE) %>% 
  as.data.frame() %>% tibble::as_tibble() %>% 
  mutate(f = himmelblau(x, y)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  ggplot(mapping = aes(x = x, y = log10_f)) +
  geom_line(size = 1.25) +
  theme_bw()

### look at 3 slices, two at the two lower modes and 1 slice in between them
df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  mutate(log_f = log10(f + 1)) %>% 
  arrange(log_f) %>% 
  slice(1:4)

expand.grid(x = seq(-5, 5, length.out = 301),
            y = c(-3.27, -2.5, -1.87),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE) %>% 
  as.data.frame() %>% tibble::as_tibble() %>% 
  mutate(f = himmelblau(x, y)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  ggplot(mapping = aes(x = x, y = log10_f)) +
  geom_line(size = 1.25,
            mapping = aes(color = as.factor(y),
                          group = y)) +
  scale_color_viridis_d('y') +
  theme_bw()

### pick one of the slices and use that as the 1D himmelblau
himmelblau_1d <- function(x, yfix)
{
  himmelblau(x, y = yfix)
}

### the 1D himmelblau function to use

yfix_use <- -2.5

tibble::tibble(
  x = seq(-5, 5, length.out = 501)
) %>% 
  mutate(f = himmelblau_1d(x, yfix = yfix_use)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  ggplot(mapping = aes(x = x, y = log10_f)) +
  geom_line(size = 1.25) +
  theme_bw()

### define the candidate grid for the horizontal coordinates
x_candidate <- tibble::tibble(
  x = seq(-5, 5, length.out = 2501)
)

### define the initial design - using evenly spaced points
x_design <- tibble::tibble(
  x = seq(-5, 5, length.out = 4)
) %>% 
  mutate(f = himmelblau_1d(x, yfix = yfix_use)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  select(x, response = log10_f)

x_design

### compare the initial design to the true function
tibble::tibble(
  x = seq(-5, 5, length.out = 501)
) %>% 
  mutate(f = himmelblau_1d(x, yfix = yfix_use)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  ggplot(mapping = aes(x = x, y = log10_f)) +
  geom_line(size = 1.25) +
  geom_point(data = x_design,
             mapping = aes(x = x, y = response),
             color = 'red', size = 5, shape = 0) +
  theme_bw()

### setup the functions we will use
library(GPfit)

### define a function which calculates the expected improvement acquisition
### function, given GP predictions
calculate_ei <- function(pred_mean, pred_sd, y_best)
{
  if( pred_sd == 0 ){
    aei <- 0
  } else {
    z_score <- (y_best - pred_mean) / pred_sd
    cdf_z <- pnorm(z_score)
    aei <- pred_sd * (z_score * cdf_z + dnorm(z_score))
  }
  
  aei
}

### define a function for range scaling the input to between 0 and 1
### because GP_fit() requires the input to be between 0 and 1
transform_input <- function(x, x_lwr, x_upr)
{
  (x - x_lwr) / (x_upr - x_lwr)
}

### set the input information to scale the inputs
input_info <- list(
  x = list(lwr = min(x_candidate$x), upr = max(x_candidate$x))
)

### define a function which reads in a design, fits a GP, predicts the search
### grid, calculates the acquisition function over the search grid
run_ei <- function(df, x_search, gp_corr, x_info)
{
  # the current best output value
  ystar <- min(df$response)
  
  # range scale the input
  df <- df %>% 
    mutate(across(.cols = 'x',
                  .fns = transform_input,
                  x_lwr = x_info$x$lwr,
                  x_upr = x_info$x$upr))
  
  x_test <- x_search %>% 
    mutate(across(.cols = 'x',
                  .fns = transform_input,
                  x_lwr = x_info$x$lwr,
                  x_upr = x_info$x$upr))
  
  # fit the GP
  gp_fit <- GP_fit(X = df$x, Y = df$response,
                   corr = gp_corr)
  
  # predict the search grid
  gp_pred <- predict.GP(gp_fit, xnew = x_test)
  
  # evaluate the acquisition function over the search grid
  ei <- purrr::map2_dbl(gp_pred$Y_hat, sqrt(gp_pred$MSE),
                        calculate_ei,
                        y_best = ystar)
  
  # package results together
  list(gp_pred_mean = gp_pred$Y_hat,
       gp_pred_sd = sqrt(gp_pred$MSE),
       ei = ei)
}

### define the information to be used for the GP kernel, use the
### SE kernel
gp_secov <- list(type = 'exponential', power = 2)

### visualize the expected improvment over the search grid given the
### initial design
x_candidate %>% 
  mutate(ei = (run_ei(x_design, x_candidate, gp_secov, input_info))$ei) %>% 
  ggplot(mapping = aes(x = x, y = ei)) +
  geom_line(size = 1.2) +
  theme_bw()

### setup functions to run multiple iterations
run_himmelbrau_1d <- function(x, yfix)
{
  log10( himmelblau_1d(x, yfix) + 1 )
}

input_info$yfix <- yfix_use
input_info$my_func <- run_himmelbrau_1d

### this was just a single iteration, so define a wrapper function which
### runs the acqusition function, finds the input value to add, generates
### the function output for that input
find_and_add_new_point <- function(df, x_search, gp_corr, x_info)
{
  # run the acqusition function
  acc_res <- x_search %>% 
    bind_cols(run_ei(df, x_search, gp_corr, x_info) %>% 
                as.data.frame() %>% tibble::as_tibble())
  
  # find the input associated with the maximum expected improvement
  new_point <- acc_res %>% 
    filter(ei == max(ei)) %>% 
    select(-ei, -gp_pred_mean, -gp_pred_sd)
  
  # run the function to generate the output at the chosen new 
  # point
  new_df <- new_point %>% mutate(response = x_info$my_func(x, x_info$yfix))
  
  # compile all results
  list(acquisition_df = acc_res,
       old_design = df,
       new_design = new_df)
}

### identify the selected new data point based on the initial design
one_iteration <- find_and_add_new_point(x_design, x_candidate,
                                        gp_secov, input_info)

one_iteration

### compare the identified point to the true function
viz_grid <- tibble::tibble(
  x = seq(-5, 5, length.out = 501)
) %>% 
  mutate(f = himmelblau_1d(x, yfix = yfix_use)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  select(x, response = log10_f)

one_iteration %>% pluck('acquisition_df') %>% 
  ggplot(mapping = aes(x = x)) +
  geom_ribbon(mapping = aes(ymin = gp_pred_mean - 2*gp_pred_sd,
                            ymax = gp_pred_mean + 2*gp_pred_sd),
              fill = 'grey') +
  geom_line(mapping = aes(y = gp_pred_mean),
            size = 1.2) +
  geom_line(data = viz_grid,
            mapping = aes(y = response),
            color = 'red', size = 1.25, linetype = 'dashed') +
  geom_point(data = x_design,
             mapping = aes(y = response),
             size = 5.5, color = 'blue') +
  geom_point(data = one_iteration %>% 
               pluck('new_design'),
             mapping = aes(y = response),
             color = 'black', size = 5.5, shape = 15) +
  geom_point(data = one_iteration %>% 
               pluck('new_design'),
             mapping = aes(y = response),
             color = 'yellow', size = 2.5) +
  theme_bw()

### define a manager function which runs the sequential design process
### adding new points for a defined number of points from the initial design
run_bayesopt <- function(df, x_search, num_iter, gp_corr, x_info)
{
  # since this is a recursive process just use a for-loop
  res <- vector(mode = 'list', length = num_iter)
  
  for(n in 1:num_iter){
    res[[n]] <- find_and_add_new_point(df, x_search, gp_corr, x_info)
    
    # add data to the design and repeat
    df <- df %>% bind_rows(res[[n]]$new_design)
  }
  
  list(all_results = res,
       final_data = df)
}

### execute the bayesian optimization for a specified number of iterations
results_from_ei <- run_bayesopt(x_design, x_candidate,
                                num_iter = 12,
                                gp_corr = gp_secov, 
                                x_info = input_info)

results_from_ei %>% glimpse()

### visualize the acquisition function at each iteration
purrr::map2_dfr(results_from_ei$all_results,
                seq_along(results_from_ei$all_results),
                function(ll, lid){
                  ll$acquisition_df %>% mutate(iter = lid)
                }) %>% 
  ggplot(mapping = aes(x = x, y = ei)) +
  geom_line(size = 1.) +
  geom_vline(data = purrr::map2_dfr(results_from_ei$all_results,
                                    seq_along(results_from_ei$all_results),
                                    function(ll, lid){
                                      ll$new_design %>% mutate(iter = lid)
                                    }),
             mapping = aes(xintercept = x),
             color = 'red', linetype = 'dashed', size = 1.) +
  facet_wrap(~iter, labeller = 'label_both') +
  theme_bw()

### visualize the GP predictions and the data ponts added to the 
### training set by the Bayesian Optimization algorithm after each iteration
purrr::map2_dfr(results_from_ei$all_results,
                seq_along(results_from_ei$all_results),
                function(ll, lid){
                  ll$acquisition_df %>% mutate(iter = lid)
                }) %>% 
  ggplot(mapping = aes(x = x)) +
  geom_ribbon(mapping = aes(ymin = gp_pred_mean - 2*gp_pred_sd,
                            ymax = gp_pred_mean + 2*gp_pred_sd,
                            group = iter),
              fill = 'grey') +
  geom_line(mapping = aes(y = gp_pred_mean,
                          group = iter),
            size = 1.2) +
  geom_point(data = purrr::map2_dfr(results_from_ei$all_results,
                                    seq_along(results_from_ei$all_results),
                                    function(ll, lid){
                                      ll$old_design %>% mutate(iter = lid)
                                    }),
             mapping = aes(x = x, y = response),
             color = 'blue', size = 3) +
  geom_point(data = purrr::map2_dfr(results_from_ei$all_results,
                                    seq_along(results_from_ei$all_results),
                                    function(ll, lid){
                                      ll$new_design %>% mutate(iter = lid)
                                    }),
             mapping = aes(y = response),
             color = 'black', size = 3, shape = 15) +
  geom_line(data = viz_grid,
            mapping = aes(y = response),
            color = 'red', size = 1, linetype = 'dotted') +
  geom_point(data = purrr::map2_dfr(results_from_ei$all_results,
                                    seq_along(results_from_ei$all_results),
                                    function(ll, lid){
                                      ll$new_design %>% mutate(iter = lid)
                                    }),
             mapping = aes(y = response),
             color = 'yellow', size = 1.1) +
  facet_wrap(~iter, labeller = 'label_both') +
  labs(y = 'y') +
  theme_bw()

### look at the expected improvement associated with the identified new
### point across iterations
purrr::map2_dfr(results_from_ei$all_results,
                seq_along(results_from_ei$all_results),
                function(ll, lid){
                  ll$acquisition_df %>% mutate(iter = lid)
                }) %>% 
  group_by(iter) %>% 
  mutate(max_ei = max(ei)) %>% 
  ungroup() %>% 
  filter(ei == max_ei) %>% 
  ggplot(mapping = aes(x = iter, y = max_ei)) +
  geom_line(size = 1.25) +
  geom_point(size = 4.5) +
  labs(x = 'iteration', y = 'Maximum EI per iteration') +
  theme_bw()

### look at the identified global minimum response compared to the 
### observed output at each iteration
results_from_ei$all_results[[1]]$old_design %>% 
  mutate(iter = 0) %>% 
  bind_rows(purrr::map2_dfr(results_from_ei$all_results,
                            seq_along(results_from_ei$all_results),
                            function(ll, lid){
                              ll$new_design %>% mutate(iter = lid)
                            })) %>% 
  mutate(global_min = cummin(response)) %>% 
  ggplot(mapping = aes(x = iter, y = response)) +
  stat_summary(geom = 'linerange',
               fun.min = 'min', fun.max = 'max',
               mapping = aes(group = iter),
               color = 'grey50', size = 1.25) +
  geom_line(mapping = aes(x = iter, y = global_min),
            color = 'blue', size = 1.25) +
  geom_point(size = 4) +
  # coord_cartesian(ylim = c(0, max(results_from_ei$final_data$response)+0.085)) +
  theme_bw()
