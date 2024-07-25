Directory containing the analysis of the results from the parameters definition and the computation of the specularity maps.

### 1_Simulations
Results analysis on the simulation performed on SIMSONIC
##GenerateMap
FOR GENERATED MAP 
- ComputeResults : programme to compute the Map, SpecularMap and metrics from the postProcess data.
- DisplayResultsSimulations : exploratory programme to find ways to display the results.
- getMetricsROI : function that return the metrics from specularity maps.
- ResultsAll : matlab files containing the Results from all the simulations performed during the internship on generated map (see calculSSD/salome).
	
##ExVivoMap
FOR EX-VIVO MAP
- ComputeResults : programme to compute the Map, SpecularMap and metrics from the postProcess data in the simulations from ex-vivo images.
- metricSimulationExVivo : exploratory programme that return the metrics from specularity maps, and define the ROI.
- ResultsAll : matlab files containing the Results from all the simulations performed during the internship one ex-vivo map (see calculSSD/salome).

###2_ExVivo
Result analysis on the experimentation performed with the Verasonic by Amadou S. Dia.
- getFig : this programme extract the region of interest on the specular images and compute the specular metrics.
- getParameters : this programme computes the parameter for each zone of the bone sample, on one image every 9 micrometers. 
- getParametersAvg : this programme computes the parameter for each zone of the bone sample, on the averaged images. 
- getCorrelationPearson : this programme compute the pearson coefficient between the parameters and the specular metrics.
- displayAveraged : this programme display the specular map and the averaged image of the bone on all zone.
##Data
- 227G, 245D, 267G, 271G : contains the averaged image in each zone, the specular map, the segmented endosteum and the beamformed images
- boneParameters : matlab file containing the microstucture parameters for each zone of the bone.
- boneParametersKc : matlab file containing the microstucture parameters for each zone of the bone, the cutt off frequency of the filter that allows to compute the surface parameters is computed for each slice, thus more accurate. Contains additionnaly the 271G bone, more degraded.
- boneSpeed : matlab file containing the speed of sound in each zone of the bone, from Amadou S. Dia.
- experimentalData : matlab file containing the specular metrics for each zone of the bone, from Amadou S. Dia.

###3_Report
- Results : live script to display image and initial work.
- Figure : COntains all the figures used in the poster/presentation/report.
