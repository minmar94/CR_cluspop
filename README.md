# Capture-recapture models for clustered populations via finite-mixtures

 1) **BUGS models** contains the following four JS-type models written in BUGS language:

     a) _JSvanilla.txt_: classical JS model (Royle and Dorazio 2020) 
      with homogeneous detection and survival probabilities;
   
     b) _JSpmix.txt_: Pledger et al. (2003)'s model with G class-specific detection probabilities, i.e. $p_1,\dots,p_G$;

     c) _JSphipmix.txt_: Pledger et al. (2003)'s model with G class-specific couple of parameters, i.e. ($p_1,\phi_1$),$\dots$,($p_G,\phi_G$);

     d) _RPT.txt_: RPT parsimonious model for characterising a population such as the one described in the 'Motivation' of the paper.

 2) **simulations** containing an $\texttt{R}$ script for carrying out the simulation study of settings A and B. The two RData files contain tables which summarises the main characteristics of each scenario.

 3) **realdata** containing the $\texttt{R}$ script for preparing and fitting JS-type models to Bottlenose dolphins' real data ("dolphin_4years.csv").

 4) **fun.R**: $\texttt{R}$ file containing functions to simulate and fit JS-type models
