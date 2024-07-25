Directory containing the toolbox and post processing programmes to run a simulation of bone insonifed by ultrasound in order to reconstruct an image using a specular transform.

### ProgrammeMatlab
Contains the useful programme to generate the necessary files to fo a SIMSONIC simulation and to post-process them in order to obtain the specular maps.

##GenerateSimulationFiles
#GenerateMap 
- MakeFile : programme used to generate all the files mandatory for the simulation in the calculSSD directory.
GenerateParameter : function to modify to generate the wanted parameters for the simulation.
- Multiples : Same but to generate all at once

- MakeGeometry : Contains the programme to generate the 2D map. 
	**MakeGeometry : generates a map representing the bone/soft tissue interface based on the provided parameters. (USEFUL)
	MakeGeometryRealistic : Generate realistic map regarding the type of bone (young, aged, osteoporotic).
	PlotFilter : filttration of the map generated to even more realisme
	
- MakeSgl : To generate the impulsed signal. 

#ExVivoMap
- **MakeGeometryExVivo : generate map from CT-images. (USEFUL)
- MakeGeometrySegmented : generate map with interface as the segmented profile.

Same organisation for generated parameters.

#GeometryMap
This files contains the geometrymap of the generated map, as it is random you can not recreate the replica therefore the original one are stored here.

##PostProcessing 
Contains the function and programme used to post-process the output of the SIMSONIC simulations.

- getPostProcess : function to compute all the post-process files, using the functions in the Specularity directory.

- LoadRfData : gather data in a recorded parameter in order to compute the specularity. 

- mainPostProcessAll : allow to post-process all simulations at once.
- mainPostProcess : post-process one simulation, with details. 

- VisualizeSimulation : programme to plot the snapshot film.

- Specularity -- directory containing the function step by step to comput the specularity map and the sub-function from Amadou S. Dia. 

##RunSimulation
Bash commande or Matlab programm to launch the SIMSONIC simulations
- run_simul : To launch the 96 simulation for each elements 
- run_simul_modify : To launch certains simulation for elements that didn't work or else.
- RunAllSimu : Matlab programme to launch simulations
- BashCommande : For copy/past commande in the Terminal 
CopiePostProcess : in order to copie a simulation from calculSSD
