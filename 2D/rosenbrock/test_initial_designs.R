### 2D example with the Rosenbrock function

library(tidyverse)

### define the function just for 2 variables
rosenbrock <- function(x, y, a = 1, b = 100)
{
  (a - x)^2 + b*(y - x^2)^2
}

### create a visualization grid
df_viz <- expand.grid(x1 = seq(-3, 3, length.out = 251),
                      x2 = seq(-2, 4, length.out = 251),
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE) %>% 
  as.data.frame() %>% tibble::as_tibble()

### original scale
df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_raster(mapping = aes(fill = f)) +
  geom_contour(mapping = aes(z = f),
               color = 'white') +
  scale_fill_viridis_c() +
  theme_bw()

df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  summary()

df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  ggplot(mapping = aes(x = f)) +
  geom_histogram(bins = 51) +
  theme_bw()

### log10 scale
df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_raster(mapping = aes(fill = log10(f))) +
  geom_contour(mapping = aes(z = log10(f)),
               color = 'white') +
  scale_fill_viridis_c() +
  theme_bw()

df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  ggplot(mapping = aes(x = log10(f))) +
  geom_histogram(bins = 51) +
  theme_bw()

### identify the locations closest to the global minimum
df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  mutate(f_log = log10(f)) %>% 
  mutate(f_min = min(f_log)) %>% 
  mutate(near_minimum = f_log <= -1) %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_raster(mapping = aes(fill = near_minimum)) +
  # geom_contour(mapping = aes(z = log10(f)),
  #              color = 'white') +
  coord_equal() +
  scale_fill_brewer("Near global minimum", palette = "Set1") +
  theme_bw()

### we need an initial design before we can execute the Bayesian optimization

### first try a simple full-factorial grid
grid_design_A <- expand.grid(x1 = seq(-3, 3, length.out = 5),
                             x2 = seq(-2, 4, length.out = 5),
                             KEEP.OUT.ATTRS = FALSE,
                             stringsAsFactors = FALSE) %>% 
  as.data.frame() %>% tibble::as_tibble() %>% 
  mutate(f = rosenbrock(x1, x2))

### overlay the initial design on the complete visualization picture
df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  mutate(f_log = log10(f)) %>% 
  mutate(f_min = min(f_log)) %>% 
  mutate(near_minimum = f_log <= -1) %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_raster(mapping = aes(fill = near_minimum)) +
  geom_point(data = grid_design_A,
             shape = 0, color = 'black', size = 7) +
  coord_equal() +
  scale_fill_brewer("Near global minimum", palette = "Set1") +
  theme_bw()

### look at the response vs the inputs in the initial design
grid_design_A %>% 
  ggplot(mapping = aes(x = x1, y = log10(f))) +
  geom_point(size = 7,
             mapping = aes(color = as.factor(x2))) +
  geom_line(size = 1.2,
            mapping = aes(color = as.factor(x2))) +
  scale_color_viridis_d("x2") +
  theme_bw()

### what if we had a smaller grid, 3x3?
grid_design_B <- expand.grid(x1 = seq(-3, 3, length.out = 3),
                             x2 = seq(-2, 4, length.out = 3),
                             KEEP.OUT.ATTRS = FALSE,
                             stringsAsFactors = FALSE) %>% 
  as.data.frame() %>% tibble::as_tibble() %>% 
  mutate(f = rosenbrock(x1, x2))

df_viz %>% 
  mutate(f = rosenbrock(x1, x2)) %>% 
  mutate(f_log = log10(f)) %>% 
  mutate(f_min = min(f_log)) %>% 
  mutate(near_minimum = f_log <= -1) %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_raster(mapping = aes(fill = near_minimum)) +
  geom_point(data = grid_design_B,
             shape = 0, color = 'black', size = 7) +
  coord_equal() +
  scale_fill_brewer("Near global minimum", palette = "Set1") +
  theme_bw()

grid_design_B %>% 
  ggplot(mapping = aes(x = x1, y = log10(f))) +
  geom_point(size = 7,
             mapping = aes(color = as.factor(x2))) +
  geom_line(size = 1.2,
            mapping = aes(color = as.factor(x2))) +
  scale_color_viridis_d("x2") +
  theme_bw()

### space filling designs are NOT full-factorial create LHS designs
library(lhs)

set.seed(222)
lhs_rand_A <- randomLHS(n = 5, k = 2) %>% 
  as.matrix() %>% as.data.frame() %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("x1", "x2"))

### look at the LHS design in the [0, 1] scale
lhs_rand_A %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_point(size = 7) +
  geom_vline(xintercept = seq(0, 1, length.out = 6), color='grey', size =1.2) +
  geom_hline(yintercept = seq(0, 1, length.out = 6), color='grey', size =1.2) +
  coord_equal() +
  theme_bw()

### random design with more points
set.seed(555)
lhs_rand_B <- randomLHS(n = 9, k = 2) %>% 
  as.matrix() %>% as.data.frame() %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("x1", "x2"))

lhs_rand_B %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_point(size = 7) +
  geom_vline(xintercept = seq(0, 1, length.out = 10), color='grey', size =1.2) +
  geom_hline(yintercept = seq(0, 1, length.out = 10), color='grey', size =1.2) +
  coord_equal() +
  theme_bw()

### iinstead of a random LHS, we can maximize the distance between
### all points while preserving the Latin property
set.seed(222)
lhs_maxmin_A <- maximinLHS(n = 5, k = 2, maxIter = 5001, eps=0.01,
                           optimize.on = 'result', method = 'iterative') %>% 
  as.matrix() %>% as.data.frame() %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("x1", "x2"))

lhs_maxmin_A %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_point(size = 7) +
  geom_vline(xintercept = seq(0, 1, length.out = 6), color='grey', size =1.2) +
  geom_hline(yintercept = seq(0, 1, length.out = 6), color='grey', size =1.2) +
  coord_equal() +
  theme_bw()


### try again
set.seed(333)
lhs_maxmin_B <- maximinLHS(n = 5, k = 2, maxIter = 5001, eps=0.01,
                           optimize.on = 'result', method = 'iterative') %>% 
  as.matrix() %>% as.data.frame() %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("x1", "x2"))

lhs_maxmin_B %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_point(size = 7) +
  geom_vline(xintercept = seq(0, 1, length.out = 6), color='grey', size =1.2) +
  geom_hline(yintercept = seq(0, 1, length.out = 6), color='grey', size =1.2) +
  coord_equal() +
  theme_bw()

### try a larger number of initial design points
set.seed(444)
lhs_maxmin_C <- maximinLHS(n = 9, k = 2, maxIter = 5001, eps=0.01,
                           optimize.on = 'result', method = 'iterative') %>% 
  as.matrix() %>% as.data.frame() %>% 
  tibble::as_tibble() %>% 
  purrr::set_names(c("x1", "x2"))

lhs_maxmin_C %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_point(size = 7) +
  geom_vline(xintercept = seq(0, 1, length.out = 10), color='grey', size =1.2) +
  geom_hline(yintercept = seq(0, 1, length.out = 10), color='grey', size =1.2) +
  coord_equal() +
  theme_bw()

### rescale from the [0,1] to the original bounds on the inputs
lhs_rand_A %>% 
  mutate(across(.cols = 'x1',
                .fns = function(x, lwr, upr){
                  (upr - lwr)*x + lwr
                },
                lwr = -3, upr = 3)) %>% 
  mutate(across(.cols = 'x2',
                .fns = function(x, lwr, upr){
                  (upr - lwr)*x + lwr
                },
                lwr = -2, upr = 4)) %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_point(size = 5, shape = 16) +
  geom_point(data = grid_design_A,
             size = 7, shape = 0) +
  coord_equal() +
  theme_bw()

### a different LHS design
lhs_maxmin_C %>% 
  mutate(across(.cols = 'x1',
                .fns = function(x, lwr, upr){
                  (upr - lwr)*x + lwr
                },
                lwr = -3, upr = 3)) %>% 
  mutate(across(.cols = 'x2',
                .fns = function(x, lwr, upr){
                  (upr - lwr)*x + lwr
                },
                lwr = -2, upr = 4)) %>% 
  ggplot(mapping = aes(x = x1, y = x2)) +
  geom_point(size = 5, shape = 16) +
  geom_point(data = grid_design_A,
             size = 7, shape = 0) +
  coord_equal() +
  theme_bw()
