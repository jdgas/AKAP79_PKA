# Aproximate Bayesian Computing (ABC) for parameter estimation in the AKAP79 -PKA model

This repository contains source code and supplementary figures to the publication
```
Church, Timothy W., Parul Tewatia, Saad Hannan, João Antunes, Olivia Eriksson, 
Trevor G. Smart, Jeanette Hellgren Kotaleski, and Matthew Gold. 
"AKAP79 enables calcineurin to directly suppress protein kinase A activity." bioRxiv (2021).
```
The ABC methodology used was based on the the study:
```
Eriksson, Olivia, Alexandra Jauhiainen, Sara Maad Sasane, Andrei Kramer, Anu G. Nair, 
Carolina Sartorius, and Jeanette Hellgren Kotaleski. "Uncertainty quantification, 
propagation and characterization by Bayesian analysis combined with global sensitivity analysis 
applied to dynamical intracellular pathway models." Bioinformatics 35, no. 2 (2019): 284-292.
````

and the corresponding code repository:
https://github.com/alexjau/uqsa

This repository is under construction and will contain the following subfolders:
- **Experimental data**, the experimental wilde type data used for model calibration as well as the experimental data from the mutation experiments.
- **Models**, the small and extended model in Matlabs Simbiology format. Five different parameter sets that fit data will be included. 
- **Posterior distribution**, the sampled parameter set fitting the experimental data.
- **Supplementary figures**, i) modelling of responses at different cAMP levels ii) illustration of five different parameter sets.
- **UQ-AKAP79**, the modified ABC code for the AKAP79-PKA model in R.
- **Figure reproduction**, code for reproducing figures in matlab.
