### optimize the 2D himmelblau function using bayesian optimization

### use a 3x3 full factorial design as the initial design

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

df_viz %>% dim()

df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = f)) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  summary()

### there is one optimum with x=3 and y=2
df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  filter(f == min(f))

### however there are 3 others with the exact same function value of 0
### these other 3 optimums are not integers
# x=-2.805118, y=3.131312
# x=-3.779310, y=-3.283186
# x=3.584428, y=-1.848126

### the visualization grid probably gets close but not exact to the 3 
### other optimums
df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  filter(abs(f - min(f)) < 3e-2)

### visualize the output in the log-scale
df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = log10(f + 1))) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

### find the 4 smallest values in the log scale
df_viz %>% 
  mutate(f = himmelblau(x, y)) %>% 
  mutate(log_f = log10(f + 1)) %>% 
  arrange(log_f) %>% 
  slice(1:4)

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

### test out a simple full-factorial grid as the initial design
init_grid <- expand.grid(x = seq(-4.75, 4.75, length.out = 3),
                         y = seq(-4.75, 4.75, length.out = 3),
                         KEEP.OUT.ATTRS = FALSE,
                         stringsAsFactors = FALSE) %>% 
  as.data.frame() %>% tibble::as_tibble()

### overlay the initial design on the surface
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
  geom_point(data = init_grid,
             size = 5, shape = 3, color = 'red') +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

### range scale the inputs to 0 and 1 to support GPfit
### recipes with step_range() unfortunately truncates values outside
### the training bounds to 0 and 1 instead of range scaling to
### negative numbers and values greater than 1!!!
### so must do the ranging scaling manually

### calculate the input min and max on the design
init_lwr <- purrr::map_dbl(init_grid, min)

init_upr <- purrr::map_dbl(init_grid, max)

init_lwr

init_upr

### range scale the initial design
norm_grid <- purrr::pmap_dfc(list(init_grid, init_lwr, init_upr),
                             function(x, xl, xu)
                             {
                               (x - xl) / (xu - xl)
                             })

### range scale the visualization grid
norm_viz <- purrr::pmap_dfc(list(df_viz, init_lwr, init_upr),
                            function(x, xl, xu){
                              (x - xl) / (xu - xl)
                            })

norm_viz %>% summary()

### plot the initial design in the ranged scaled inputs and the true function
norm_viz %>% 
  bind_cols(df_viz %>% mutate(f = himmelblau(x, y)) %>% 
              select(f)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = log10_f)) +
  geom_point(data = norm_grid,
             color = 'red', shape = 3, size = 5) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

### before running the full bayesian optimization, execute a single iteration

### the starting design with the log-transformed output
df <- norm_grid %>% 
  bind_cols(init_grid %>% 
              mutate(f = himmelblau(x, y)) %>% 
              mutate(response = log10(f + 1)) %>% 
              select(response))

df

### start by specifying the GP information

library(GPfit)

### define the information to be used for the GP kernel, use the
### SE kernel
gp_secov <- list(type = 'exponential', power = 2)

# fit the GP using the normalized input grid
set.seed(111)
gp_fit <- GP_fit(X = df %>% select(-response) %>% as.matrix(), 
                 Y = df$response,
                 corr = gp_secov,
                 maxit = 1001)

### the GP_fit() result seems like it can be somewhat "unstable"
### may not get the same result everytime we run the algorithm
### might need to switch to dicekriging

### to try and make the results more consistent I set the seed
### and raised the maxit to 1001 from the defeault of 100

# predict the search grid
gp_pred <- predict.GP(gp_fit, xnew = norm_viz %>% as.matrix())

### visualize the predicted mean and predicted sd as a surface
### on the normalized input visualization grid
norm_viz %>% 
  mutate(pred_mean = gp_pred$Y_hat) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = pred_mean)) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

norm_viz %>% 
  mutate(pred_sd = sqrt(gp_pred$MSE)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = pred_sd)) +
  coord_equal() +
  scale_fill_viridis_c(option = 'magma') +
  theme_bw()

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

### calculate the expected improvement for each prediction point
current_best <- min(df$response)

test_ei <- purrr::map2_dbl(gp_pred$Y_hat, sqrt(gp_pred$MSE),
                           calculate_ei,
                           y_best = current_best)

### visualize the expected improvement as a surface
norm_viz %>% 
  mutate(EI = test_ei) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = EI)) +
  coord_equal() +
  scale_fill_viridis_c(option = 'plasma') +
  theme_bw()

### the new point is selected as the point with the highest expected 
### improvement
new_point <- norm_viz %>% 
  mutate(EI = test_ei) %>% 
  filter(EI == max(EI))

### add the point on top of the initial design
norm_viz %>% 
  bind_cols(df_viz %>% mutate(f = himmelblau(x, y)) %>% 
              select(f)) %>% 
  mutate(log10_f = log10(f + 1)) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_raster(mapping = aes(fill = log10_f)) +
  geom_point(data = norm_grid,
             color = 'red', shape = 3, size = 5) +
  geom_point(data = new_point,
             color = 'magenta', shape = 4, size = 5) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_bw()

### define a function to transform the inputs
transform_inputs <- function(x, xlwr, xupr)
{
  (x - xlwr) / (xupr - xlwr)
}

### define a function to undo the transformation
backtransform_inputs <- function(u, xlwr, xupr)
{
  (xupr - xlwr) * u + xlwr
}

### define a function which reads in a design, fits a GP, predicts the search
### grid, calculates the acquisition function over the search grid
run_ei <- function(my_df, x_search, gp_corr, my_info)
{
  x_design <- my_df %>% select(all_of(my_info$input_names))
  
  response <- my_df %>% select(all_of(my_info$response_name)) %>% pull()
  
  # current best output value
  ystar <- min(response)
  
  # range scale the inputs
  x_ready <- purrr::pmap_dfc(list(x_design, my_info$x_lwr, my_info$x_upr),
                             my_info$transform_inputs)
  
  # in case the search grid inputs are out of order
  x_search <- x_search %>% select(all_of(my_info$input_names))
  
  x_test <- purrr::pmap_dfc(list(x_search, 
                                 my_info$x_lwr, 
                                 my_info$x_upr),
                            my_info$transform_inputs)
  
  # standardize current outputs
  y_ready <- (response - mean(response)) / sd(response)
  # y_ready <- response
  
  # browser()
  
  # fit the GP
  set.seed(my_info$seed_use)
  gp_fit <- GP_fit(X = as.matrix(x_ready), 
                   Y = y_ready,
                   corr = gp_corr,
                   maxit = my_info$maxit)
  
  # predict the search grid
  gp_pred <- predict.GP(gp_fit, xnew = as.matrix(x_test))
  
  # undo the output standardization
  pred_mean <- sd(response) * gp_pred$Y_hat + mean(response)
  pred_sd <- sd(response) * sqrt(gp_pred$MSE)
  # pred_mean <- gp_pred$Y_hat
  # pred_sd <- sqrt(gp_pred$MSE)
  
  # evaluate the acquisition function over the search grid
  ei <- purrr::map2_dbl(pred_mean, pred_sd,
                        calculate_ei,
                        y_best = ystar)
  
  # package results together
  list(gp_pred_mean = pred_mean,
       gp_pred_sd = pred_sd,
       ei = ei)
}

### wrap the himmelbrau function so that it returns the log-transformed value
run_himmelbrau <- function(x, y)
{
  log10( himmelblau(x, y) + 1 )
}

### this was just a single iteration, so define a wrapper function which
### runs the acqusition function, finds the input value to add, generates
### the function output for that input
find_and_add_new_point <- function(my_design, x_search, gp_corr, my_info)
{
  # run the acqusition function
  acc_res <- x_search %>% 
    bind_cols(run_ei(my_design, x_search, gp_corr, my_info) %>% 
                as.data.frame() %>% tibble::as_tibble())
  
  # find the input associated with the maximum expected improvement
  # if multiple points are associated with the same max acquisition
  # function, for simplicity just try the first 11
  new_point_results <- acc_res %>% filter(ei == max(ei)) %>% 
    tibble::rowid_to_column('rowid_id') %>% 
    filter(rowid_id < 12) %>% 
    select(-rowid_id)
  
  new_point <- new_point_results %>% 
    select(all_of(my_info$input_names))
  
  # run the function to generate the output at the chosen new 
  # point
  new_df <- new_point %>% 
    rowwise() %>% 
    mutate(response = my_info$the_func(x, y)) %>% 
    ungroup()
  
  # compile all results
  list(acquisition_df = new_point_results,
       old_design = my_design,
       new_design = new_df)
}

### compile the necessary info together
info_use_01 <- list(
  x_lwr = init_lwr,
  x_upr = init_upr,
  the_func = run_himmelbrau,
  transform_inputs = transform_inputs,
  backtransform_inputs = backtransform_inputs,
  input_names = c("x", 'y'),
  response_name = 'response',
  maxit = 1001,
  the_seed = 56781
)

### the starting data
init_data <- init_grid %>% mutate(response = run_himmelbrau(x, y))

### find 1 point
one_iter_res <- find_and_add_new_point(init_data, df_viz, gp_secov, info_use_01)

### gp prediction and EI for that one result
one_iter_res$acquisition_df

### define a manager function which runs the sequential design process
### adding new points for a defined number of points from the initial design
run_bayesopt <- function(my_df, x_search, num_iter, gp_corr, x_info)
{
  # since this is a recursive process just use a for-loop
  res <- vector(mode = 'list', length = num_iter)
  
  for(n in 1:num_iter){
    res[[n]] <- find_and_add_new_point(my_df, x_search, gp_corr, x_info)
    
    # add data to the design and repeat
    my_df <- my_df %>% bind_rows(res[[n]]$new_design)
  }
  
  # added_df <- purrr::map_dfr(res, 'new_design')
  
  list(all_results = res,
       final_data = my_df)
}

### execute the bayesian optimization for a specified number of iterations
results_from_ei <- run_bayesopt(init_data, df_viz,
                                num_iter = 51,
                                gp_corr = gp_secov, 
                                x_info = info_use_01)

### visualize the acquisitioin function at each iteration
results_from_ei$all_results %>% length()

results_from_ei$all_results[[1]] %>% class()

results_from_ei$all_results[[1]]$acquisition_df %>% glimpse()

### look at how the design grows at each iteration and which point is added
### include the global optimum points as reference
purrr::map2_dfr(results_from_ei$all_results,
                seq_along(results_from_ei$all_results),
                function(ll, lid){
                  ll$old_design %>% mutate(iter = lid) %>% 
                    mutate(type = 'design')
                }) %>% 
  bind_rows(purrr::map2_dfr(results_from_ei$all_results,
                            seq_along(results_from_ei$all_results),
                            function(ll, lid){
                              ll$new_design %>% mutate(iter = lid) %>% 
                                mutate(type = 'added')
                            })) %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  geom_point(mapping = aes(color = type, shape = type)) +
  geom_point(data = df_viz %>% 
               mutate(f = himmelblau(x, y)) %>% 
               mutate(log_f = log10(f + 1)) %>% 
               arrange(log_f) %>% 
               slice(1:4),
             shape = 0, color = 'red', size = 3) +
  coord_equal() +
  facet_wrap(~iter, labeller = 'label_both') +
  scale_color_manual('',
                     values = c("added" = 'blue',
                                'design' = 'grey50')) +
  scale_shape_manual('',
                     values = c("added" = 4,
                                "design" = 16)) +
  theme_bw()

### plot the output value observed at each iteration and compare to the
### currently identified global minimum
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
  coord_cartesian(ylim = c(0, max(results_from_ei$final_data$response)+0.085)) +
  theme_bw()

### plot the GP prediction and the observed response for each iteration
results_from_ei$all_results[[1]]$old_design %>% 
  mutate(iter = 0) %>% 
  bind_rows(purrr::map2_dfr(results_from_ei$all_results,
                            seq_along(results_from_ei$all_results),
                            function(ll, lid){
                              ll$new_design %>% mutate(iter = lid)
                            })) %>% 
  mutate(global_min = cummin(response)) %>% 
  left_join(purrr::map2_dfr(results_from_ei$all_results,
                            seq_along(results_from_ei$all_results),
                            function(ll, lid){
                              ll$acquisition_df %>% 
                                select(gp_pred_mean, gp_pred_sd, ei) %>% 
                                mutate(iter = lid)
                            }),
            by = c("iter")) %>% 
  ggplot(mapping = aes(x = iter)) +
  stat_summary(geom = 'linerange',
               fun.min = 'min', fun.max = 'max',
               mapping = aes(group = iter,
                             y = response),
               color = 'grey50', size = 1.25) +
  geom_line(mapping = aes(x = iter, y = global_min),
            color = 'blue', size = 1.25) +
  geom_linerange(mapping = aes(ymin = gp_pred_mean - gp_pred_sd,
                               ymax = gp_pred_mean + gp_pred_sd,
                               group = iter),
                 size = 1, color = 'darkorange', alpha = 0.55) +
  geom_point(mapping = aes(y = gp_pred_mean),
             color = 'darkorange', size = 2.5) +
  geom_point(size = 4,
             mapping = aes(y = response)) +
  coord_cartesian(ylim = c(0, max(results_from_ei$final_data$response)+0.085)) +
  theme_bw()

### plot the acquisition function with respect to the iteration
purrr::map2_dfr(results_from_ei$all_results,
                seq_along(results_from_ei$all_results),
                function(ll, lid){
                  ll$acquisition_df %>% mutate(iter = lid)
                }) %>% 
  ggplot(mapping = aes(x = iter, y = ei)) +
  geom_line(size = 1.2) +
  geom_point(size = 5) +
  theme_bw()

### 