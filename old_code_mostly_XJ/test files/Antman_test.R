library(AntMAN)

mcmc = AM_mcmc_parameters(niter=20000)
mix = AM_mix_hyperparams_multinorm ()
prior_component = AM_mix_components_prior_dirac(10) # 10 colours present
fit = AM_mcmc_fit (AM_sample_unipois()$y, prior_component, mcmc)
summary (fit)


fit = AM_mcmc_fit( AM_sample_unipois()$y,
             AM_mix_hyperparams_unipois (alpha0=2, beta0=0.2),
             mcmc_parameters = AM_mcmc_parameters(niter=50, burnin=0, thin=1, verbose=0))



data(galaxy, package = "AntMAN")
y_uvn = galaxy
mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=20.83146, k0=0.3333333, nu0=4.222222, sig02=3.661027)

mcmc_params        = AntMAN::AM_mcmc_parameters(niter=20000, burnin=5000, thin=10, verbose=1)
components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AntMAN::AM_mcmc_fit(
  y = y_uvn,
  mix_kernel_hyperparams = mixture_uvn_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)

summary (fit)
plot (fit)
