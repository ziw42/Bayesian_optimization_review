### try out Bayesian optimization with Expected Improvement acquisition
### function on a multimodal test function

### the functions are setup to support visualizing the results, not necessaily
### running and storing everything as efficiently as possible

library(tidyverse)

library(GPfit)

### test function comes from machine learning mastery
# https://machinelearningmastery.com/1d-test-functions-for-function-optimization/

### define the function we wish to minimize
my_func <- function(x){sin(x) + sin((10/3)*x)}

### visualize the function over the desired bounds of the problem
viz_grid <- tibble::tibble(
  x = seq(-2.7, 7.5, length.out = 3001)
) %>% 
  mutate(f = my_func(x))

viz_grid %>% 
  ggplot(mapping = aes(x = x, y = f)) +
  geom_line(size = 1.2) +
  geom_vline(data = viz_grid %>% 
               filter(f == min(f)),
             mapping = aes(xintercept = x),
             color = 'red', linetype = 'dashed', size = 1.) +
  theme_bw()

### create the initial design
df_design <- tibble::tibble(
  x = seq(-2.7, 7.5, length.out = 4)
) %>% 
  mutate(y = my_func(x))

### compare the initial design to the true function
viz_grid %>% 
  ggplot(mapping = aes(x = x, y = f)) +
  geom_line(size = 1.2) +
  geom_vline(data = viz_grid %>% 
               filter(f == min(f)),
             mapping = aes(xintercept = x),
             color = 'red', linetype = 'dashed', size = 1.) +
  geom_point(data = df_design,
             mapping = aes(x = x, y = y),
             color = 'blue', size = 7.5) +
  theme_bw()

### identify the search grid that the GP model will search over
### to identify new points
x_candidate <- tibble::tibble(
  x = seq(-2.7, 7.5, by = 0.01)
)

x_candidate %>% nrow()

viz_grid %>% 
  ggplot(mapping = aes(x = x, y = f)) +
  geom_line(size = 1.2, color = 'grey') +
  geom_vline(data = viz_grid %>% 
               filter(f == min(f)),
             mapping = aes(xintercept = x),
             color = 'red', linetype = 'dashed', size = 1.) +
  geom_point(data = df_design,
             mapping = aes(x = x, y = y),
             color = 'blue', size = 7.5) +
  geom_point(data = x_candidate %>% 
               mutate(f = my_func(x)),
             color = 'black', size = 2) +
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

### define a function for range scaling the input to between 0 and 1
### because GP_fit() requires the input to be between 0 and 1
transform_input <- function(x, x_lwr, x_upr)
{
  (x - x_lwr) / (x_upr - x_lwr)
}

### set the input information to scale the inputs
input_info <- list(
  x = list(lwr = min(viz_grid$x), upr = max(viz_grid$x))
)

### define a function which reads in a design, fits a GP, predicts the search
### grid, calculates the acquisition function over the search grid
run_ei <- function(df, x_search, gp_corr, x_info)
{
  # the current best output value
  ystar <- min(df$y)
  
  # range scale the input
  df <- df %>% 
    mutate(across(.cols = 'x',
                  .fns = transform_input,
                  x_lwr = x_info$x$lwr,
                  x_upr = x_info$x$upr))
  
  x_search <- x_search %>% 
    mutate(across(.cols = 'x',
                  .fns = transform_input,
                  x_lwr = x_info$x$lwr,
                  x_upr = x_info$x$upr))
  
  # fit the GP
  gp_fit <- GP_fit(X = df$x, Y = df$y,
                   corr = gp_corr)
  
  # predict the search grid
  gp_pred <- predict.GP(gp_fit, xnew = x_search)
  
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
  mutate(ei = (run_ei(df_design, x_candidate, gp_secov, input_info))$ei) %>% 
  ggplot(mapping = aes(x = x, y = ei)) +
  geom_line(size = 1.2) +
  theme_bw()

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
  new_df <- new_point %>% mutate(y = my_func(x))
  
  # compile all results
  list(acquisition_df = acc_res,
       old_design = df,
       new_design = new_df)
}

### identify the selected new data point based on the initial design
one_iteration <- find_and_add_new_point(df_design, x_candidate,
                                        gp_secov, input_info)

one_iteration %>% glimpse()

one_iteration %>% pluck('acquisition_df') %>% 
  ggplot(mapping = aes(x = x)) +
  geom_ribbon(mapping = aes(ymin = gp_pred_mean - 2*gp_pred_sd,
                            ymax = gp_pred_mean + 2*gp_pred_sd),
              fill = 'grey') +
  geom_line(mapping = aes(y = gp_pred_mean),
            size = 1.2) +
  geom_line(data = viz_grid,
            mapping = aes(y = f),
            color = 'red', size = 1.25, linetype = 'dashed') +
  geom_point(data = df_design,
             mapping = aes(y = y),
             size = 5.5, color = 'blue') +
  geom_point(data = one_iteration %>% 
               pluck('new_design'),
             mapping = aes(y = y),
             color = 'black', size = 5.5, shape = 15) +
  geom_point(data = one_iteration %>% 
               pluck('new_design'),
             mapping = aes(y = y),
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
results_from_ei <- run_bayesopt(df_design, x_candidate,
                                num_iter = 16,
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
             mapping = aes(x = x, y = y),
             color = 'blue', size = 3) +
  geom_point(data = purrr::map2_dfr(results_from_ei$all_results,
                                    seq_along(results_from_ei$all_results),
                                    function(ll, lid){
                                      ll$new_design %>% mutate(iter = lid)
                                    }),
             mapping = aes(y = y),
             color = 'black', size = 3, shape = 15) +
  geom_line(data = viz_grid,
            mapping = aes(y = f),
            color = 'red', size = 1, linetype = 'dotted') +
  geom_point(data = purrr::map2_dfr(results_from_ei$all_results,
                                    seq_along(results_from_ei$all_results),
                                    function(ll, lid){
                                      ll$new_design %>% mutate(iter = lid)
                                    }),
             mapping = aes(y = y),
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

### look at the identified global minimum output value per iteration
purrr::map2_dfr(results_from_ei$all_results,
                seq_along(results_from_ei$all_results),
                function(ll, lid){
                  ll$old_design %>% mutate(iter = lid)
                }) %>% 
  group_by(iter) %>% 
  mutate(y_min = min(y)) %>% 
  ungroup() %>% 
  filter(y == y_min) %>% 
  ggplot(mapping = aes(x = iter, y = y_min)) +
  geom_line(size = 1.25) +
  geom_point(size = 6.5) +
  labs(x = 'iteration', y = 'Identified global minimum output') +
  theme_bw()
