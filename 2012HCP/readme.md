##How to run the analysis

1. Combine all data in one hdf file and all MC in another file
2. Use `preprocessing/lumi.py` to select only data from good runs
3. `data.py` defines the locations of the hdf files and the names of the dataframe
4. `parameters.py` contains some efficiencies that need to be set
4. Validation: Run each of the scripts in the validation directory and make sure everything looks reasonable. Some of these plots and tables should probably go in the analysis documentation (todo: make the plots look nicer)
5. Make a background prediction
	
		import BackgroundFit
		BackgroundFit.do_background_fit(data, mc, ...)
	
6. These background predictions get put into the datacard manually
7. Make the “money plots” using `analysis/money_plot.py` and the results table using `analysis/make_results_table.py`. These are the fundamental physics results

## How to make datacards for SMS limits

1. All mass points for a given SMS should go into the same hdf file
2. `parameters.py` defines some selection uncertainties that need to be set
3. Run the function in `limits/signal.py` to generate a set of datacards for each SMS.