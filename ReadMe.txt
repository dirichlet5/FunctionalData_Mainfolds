Main code documentation description,

- 1_Simulation.R,
	Complete the generation of simulation data and the comparison of the effect under three different anomaly detection algorithms.

- 2_Plot.R,
	Plot the line graph of TPR and FPR based on the results of the runs in 1_Simulation.R

- 3_Plot3D.R,
	Plot the four types of anomalous data in the sphere according to their curves in 1_Simulation.R
	
- 4_FlightData.R,
	Compare the effect of two different anomaly detection algorithms in real data, trajectory data.
------------------------------------------------------------------------
Secondary code documentation description,

- data_type1.R
	Generate simulation data for exception type 1.

- data_type2.R
	Generate simulation data for exception type 2.
	
- data_type3.R
	Generate simulation data for exception type 3.
	
- data_type4.R
	Generate simulation data for exception type 4.

- ord_detect.R
	Outlier detection based on an ordinary distance detection algorithm.

- JCGS_detect.R
	Outlier detection based on an algorithm from a JCGS paper.

- outlier_detect.R
	Outlier detection based on our proposed algorithm.
	
------------------------------------------------------------------------
Folder Description,

- data,
	Contains a real track dataset, python code for plotting geospatial maps from the track data, and the dataset used in the simulation experiment.

- function,
	Contains the functions used to generate the simulation data and the functions used in the calculation of the three anomaly detection algorithms.

- result_mean,
	Mean values of TPR and FPR of the assay results after repeating the experiment with the three algorithms for the simulated data.

- result_process,
	After repeating the experiment with the simulated data under the three algorithms, the values of TPR and FPR for each calculation process of the detection results.

- result_flight,
	Detection results of the two algorithms in real trajectory data.
------------------------------------------------------------------------