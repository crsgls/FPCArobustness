# FPCA and missing data

Simulation code for the empirical simulation study from the arxiv preprint [*Functional principal component analysis as an alternative to mixed-effect models for describing sparse repeated measures in presence of missing data*]() by Ségalas C, Helmer C, Genuer R, Proust-Lima C.

-   In `simus_dropout.R`: with this file, longitudinal trajectories can be generated (with various observation size `nmes`, dropout intensity `txdo` and missing data mechanisms `missing`). On a training subset of these, linear mixed model and FPCA are estimated. Then, these are used to predict missing data in a test subset of the generated longitudinal trajectories. Longitudinal RMSEs are then computed to evaluate the ability of the model to predict the missing data, and hence the robustness to missing data, depending upon the missing data scenario. The reference RMSEs is computed from the generated complete data on which a flexible mixed model with spline trend is estimated. Note that the joint models can sometimes struggle to converge, making the code looping over attempts until convergence. To avoid this, and the long running time it implies, the *MNAR* missing data scenario can be commented out by commenting lines `380-381` and uncommenting lines `382-383`.

-   In `simus_fpca_datagen.R` and `simus_fpca.R` : with the first of these files, FPCA is estimated on longitudinal data to generate an appropriate true mean function and true principal components (stored in `dt_preds.Rdata`). Then, these are loaded and FPCA individual scores are generated to reconstruct longitudinal trajectories. Noise and missing data are added (with various observation size `nmes`, dropout intensity `txdo` and missing data mechanisms `missing`). Then, FPCA are estimated from this noisy discretized data and the estimated mean function and functional principal components can then be compared to their true counterparts, used to generate the data.

Note that these 'R' code are designed to run only one replicate of the chosen scenario. They were replicated 1000 times in the empirical simulation study.
