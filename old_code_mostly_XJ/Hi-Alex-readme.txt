Hi Alex,

This folder contains all the code for the distance project. Sorry the main things have not been written down very well. I try to explain a bit here.

There are 3 big parts: 

1. about "salso": the algorithm are similar with salso, so I call our big algorithm "salso" in the name of R files. 
-- "Salso-functions-generalformat.R": functions for our algorithm, which is salso-like framework which you can input any distance-based loss functions. 
(In "test files" folder:
  -- "salso_functions.R": functions for the original salso algorithm.
  -- "salso_W2_functions.R": a middle code version which I try to put Wasserstein distance (p=2) distance into the salso framework.)

2. about "distances": I implement the following distances so far: 
2.1 Wasserstein distance with p=2 (W2): 
-- "W2-functions.cpp": the 1d W2 distance functions.
(In "test files" folder:
  -- "mixture_normal_W2.R": preliminary version (not use it anymore).
  -- "mixture_normal_W2_location_scale.R": the test file for W2 functions.
  -- "W2-Salso-functions.R")

2.2 Pearson chi-square distance:
-- "Pearson-functions.R": exactly same as "Pearson-functions-1d.R".
-- "Pearson-functions-1d.R"
-- "Pearson-functions-2d.R"
(In "test files" folder:
  -- "Pearson-chisq-distance-test.R")

2.3 Kolmogorov–smirnov distance:
-- "KS-functions-1d.R"
-- "KS-functions-2d.R"
(In "test files" folder:
  -- "KS_test.R")

3. about "simulation": full process from data generation, Antman Bayesian inference, then use partition samples from Antman to do our algorithm.
-- (Most important) "Simulation1.R": Do N times simulations, run the algorithm based on W2, KS and Pearson distances, draw boxplot and compared with other loss functions (Binder and VI in salso paper). This file outputs very nice plots. Warning: it takes quite long time like 1-2 hours, prepare to wait when you run it.
-- "MCMC-fit-plot.R": to output some plots with MCMC average density.
-- (Important) "New-salso-KS.R": whole algorithm process with KS distance.
-- (Important) "New-salso-Pearson.R": whole algorithm process with Pearson distance.
-- (Important) "New-salso-W2.R": whole algorithm process with W2 distance.
(In "test files" folder:
  -- "Antman_test.R": try the "Antman" Rpackage.
  -- "Antman_toy.R": test my "salso" code with the salso package. They give same result.
  -- "Antman_toy_salsoW2": try W2 as a loss function preliminarily.)

That's all. If you want to have a try on these codes, you could start with "New-salso-KS.R", "New-salso-Pearson.R" and "New-salso-W2.R". These three will not take long time.

Some future directions: 
1. to develop 2d versions of KS and Pearson distances (the W2 formula is only applicable on 1d case). Actually I have a preliminary version of 2d KS and Pearson, but not sure they are 100% correct.

1.1 after you get 2d distancce, try to involve 2d distance loss function in our algorithm. It will be easy, just change the input distance function to the 2d version. But also need to check if everything goes well.

2. we could try one more distance --- Maximum Mean Discrepancy, on mixture models as a new loss function.

That's all really. Thank you Alex. Hope everything goes well. Let me know if you have any questions. I am always happy to help.
