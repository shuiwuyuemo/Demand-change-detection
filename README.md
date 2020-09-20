
#########################################################################

    These are the R codes of the paper:                                                                                
        Detecting the demand changes of bike sharing: A Bayesian hierarchical approach       
    If you wanna use these codes for research purposes, please cite this paper as well           
    Require R packages: MCMCpack, MASS, geoR, mnormt, apcluster                                   
    R version: 3.4.4                                                                                                              

#########################################################################

#----------Description of codes----------#
1. demand_change_detection_weekly.r  :  Demand change detection using the weekly-integrated bike sharing demand data
2. demand_change_detection_daily.r  :  Demand change detection using the daily-integrated bike sharing demand data
3. unify_variables.r  :  Resampling the state-specific beta and sigma2, which is used to unify the independent variables and dependent variables applied in model comparison and data segment clustering
4. model_comparison_criteria.r  :  Calculating the log marginal likelihoods, AIC, BIC and DIC of the demand change detection outcomes
5. time_segment_cluster.r  :  Clustering the data segments divided by the demand changes

#----------Description of input data----------#
1. cluster_info.csv stores the information of clusters. It has num_region rows and at least one column named "cluster_label", which stores the indexes of regions.
2. cluster_generation_final_weekly.csv stores the weekly-integrated bike sharing demands and the independent variables, which account for the changing regularity of bike sharing demands. The first column is the time index, and the next num_region columns are the demand records of the num_region regions.
3. cluster_generation_final_daily.csv stores the daily-integrated bike sharing demands and the independent variables. The description
