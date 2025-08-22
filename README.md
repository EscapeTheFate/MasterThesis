Overview of R-Files:
- Ordered-Probit-DGP.R = main file used for calculations per PC variant. One slight note: Some later simulation have overwritten some earlier one's due to same parameter configuration used and the sequencing of certain parameter does yield different results, i.e. rho = 0.3, 0.6, 0.9 yields different results than seq(from = .3, to = .9, by = .3)  
- (HomePC)-Ordered-Probit-DGP.R (variant, see above ^^^^)
- (ProfPC)-Ordered-Probit-DGP.R (variant, see above ^^^^)
- IIApropertyFigure.R = Just as the name states
- ShinyAppForPlots.R = Main file whose Shiny App was initially used for findings, but some part of the final Figures in 'Pictures' were created using the last segments
- Other_Code_Space_for_plots_and_stuff.R = Other main file for graphic creation (later segments after initialising libraries)

Overview of CSV-File Datasets:
- Simulation_Results.csv (Main or final file with results)
- DGP1_... Variants (Older Results, sometimes not entirely correct/adjusted)

Overview of Folders:
- Estimate_Collection = csv dataset files from each estimated model, including all (usually) nreps = 3000 simulation replications
- True_Variance_Estimate_Collection = csv dataset files from each single large-N (=20.000) simulation, later used for the rescaled confidence interval calculation
- Pictures = Rstudio ggplot2 Pictures intended for the Text
- Old Data = just old/outdated or other files that aren't really relevant for anything anymore, i.e. GetConvergenceInfo.R that extracts convergence info for each dataset in Estimate_Collection, but info isn't really useful 
