#if !defined(CONFIG_HPP)
#define CONFIG_HPP 1

/* regard a solution approximately optimal if it is less than (thr_frac * initial_guess) */
const double thr_frac = 0.05;

/* (2L * fd_frac) will be the stepsize used by nagd.cpp for estimating gradients */
const double fd_frac = 1e-6;

/* dictates sqrt(learning_rate) in SGD, i.e., the standard deviation of the Gaussian (Brownian) noise per time unit */
//const double noise_level = 1;

const int mu_size = 10000000;
const double dt_init = 1e-4;
const double dt_grow = 1.1;
const double dt_shrink = 0.5;
const double mu_f = 0.01;
const double h_coef = 0.05;

#endif