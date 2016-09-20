#aspi - Analysis of Symmetry of Parasitic Infections
Tools for the analysis and visualization of bilateral asymmetry in parasitic infections. 

##Example data
###Observed parasite distributions
```R
diplostomum_eyes_excl_lenses
diplostomum_lenses
```

###Simulated parasite distributions
```R
simulated_asymmetry_inconsistent_bias
simulated_left_bias_heterogeneous_proportions
simulated_left_bias_homogeneous_proportions
simulated_symmetrical_infection
```

##Examples of usage
###Replicated g-tests of goodness of fit
```R
g.test(diplostomum_eyes_excl_lenses)
```

###Exact binomial tests
```R
eb.test(diplostomum_lenses)
```

###Histogram
```R
plotHistogram(diplostomum_eyes_excl_lenses)
```

###Volcano plot
```R
plotVolcano(diplostomum_eyes_excl_lenses)
```
