# CoxAIPW 0.0.3

- CoxAIPW function now includes an additional crossFit parameter, which can be turned to FALSE to stop cross-fitting when all nuisance function estimators are root-n.

- CoxAIPW function now includes an additional augmentation parameter in replacement of the previous 'random_censoring' parameter. Augmentation with respect to either propensity score, censoring or both can be selected.

- Under non-proportional hazard model, CoxAIPW estimates a causal estimated that is a weighted average of the time-varying log hazard ratio. It also outputs x and y values for producing a smoothed plot of the time-varying log hazard ratio.

- Update reference articles in the DESCRIPTION.

---


# CoxAIPW 0.0.2

- CoxAIPW function now also outputs cumulative baseline hazard function and survival function for both groups, all of which evaluated at each follow up time.

- CoxAIPW function now supports an additional 'random_censoring' parameter, which allows simplified estimation under the random censoring assumption instead of the informative censoring assumption.

- CoxAIPW function now supports data with tied event times.

- Update referenced articles in the DESCRIPTION.

---


# CoxAIPW 0.0.1

* Added a `NEWS.md` file to track changes to the package.
