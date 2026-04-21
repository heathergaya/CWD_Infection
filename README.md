# CWD_Infection
Code for: Host-pathogen dynamics and the hidden hazards of chronic wasting disease
Authors: Heather E. Gaya, Marcelo H. Jorge, Michael J. Chamberlain, Amy V. Nalls,  Nathaniel D. Denkers, Candace K. Mathiason, Mark G. Ruder, Gino J. D’Angelo, and Richard B. Chandler

Code written by Heather Gaya 
Last modified April 21, 2026

The following files are included in this repository:
- 01_CWDSurvivalModel.R contains the code to run the analysis, written in NIMBLE
    -  TVCdat4_Nov14.rds contains all the information needed to run the model (data, constants, parameters of interest, and initial values) 
- 02_CWDVisuals.R contains code to create all figures
    - UpdatedSurvival_wSex_Jan25.rds contains parameter estimates needed for several figures  
    - LargePosterior_133.rds is needed to create Figure 3, left panels
    - LargePosterior_31.rds is needed to create Figure 3, right panels
    - Predicted_Hazard_Survival.csv contains the daily hazard (survival) for fawns in the study area
    - Prevalence_Data.csv contains all testing data from deer in the study area. This include deer with GPS collars and those with only ear tags (not used in survival analysis)
