# melange-characterization
Code for analysis of iceberg size distributions extracted from very high-resolution digital elevation models of near-terminus ice melange. After manually pre-processing DEMs (i.e., tracing termini, masking out open water in order to isolate the melange, drawing a centerline), the DEMs are partitioned into bins and iceberg size distributions are estimated from binned elevation distributions. The binned size distributions are paired with ITS_LIVE speeds extracted along the centerline to estimate buttressing generated at glacier termini by melange. The size distributions can be fit with fragmentation theory equations to assess controls on iceberg size distributions with respect to buttressing. 

The only files that you need to run, in order, to execute the code are 
1_generate_size_distrib.m: A Matlab wrapper code that will call all necessary Matlab functions in this repository to create iceberg size distributions from DEMs.
2_manually_adjust_fits.ipynb (optional!): Run 2_manually_adjust_fits.ipynb with Jupyter to use sliders to manually adjust the automated fits produced with the Matlab codes so that they are tuned to match the larger iceberg size classes in each distribution.
3_import_itslive_point_timeseries.ipynb: Run this as a Google Colab notebook to download ITS_LIVE data from regularly-spaced points along the centerline. Designed to use the centerline created during iceberg size distribution analysis in 1_generate_size_distrib.m.
4_compile_melange_characteristics.m: Pulls data from the iceberg size distributions, speed profiles, and terminus time series (delineations during melange masking plus TermPicks and/or manually-delineated terminus data for your site, if available) to create seasonal profiles of melange thickness, packing density, and speed in order to estimate melange buttressing.
5_visualize_melange_spatiotemporal_variability.m: Produces site-specific and regional plots of melange attributes. Requires customization!
6_calculate_melange_meltwater_fluxes_DEV.m: CODE IN DEVELOPMENT. Combines the melange observations to estimate freshwater fluxes from the melange.

This code was developed jointly by Ellyn Enderlin, Jukes Liu, and Michal Kopera at Boise State University. Please contact ellynenderlin@boisestate.edu or post a issue or discussion on Github if you encounter problems using any codes in this repository.
