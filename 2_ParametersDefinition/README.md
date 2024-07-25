This files contain the MATLAB programs needed to obtain the parameters from the ex-vivo images.

#1_ BoneImage
Micro-CT images. One image is chosen for each zone to compute the parameters

#2_EndosteumDefinition
This file contains function used in order to compute the parameters. The usual step are the following :

ExtractBoundary : Return the segmented endosteum as well as the pores from the binarized bone images.
FilterProfil : Filter the parabolic form of the endosteum to normalize the computation of the surface parameters. 
To compute the parameters in the near intra-cortical region of interest :
FitParabola : Return the coordinate of the parabola fitted to the endosteum.
ExpandParabola : Return the coordinate of the expended parabola of a given wifth of extention (usually one wavelength)

StudyFrequency : 
Contains programmes used in the study of the filtration of the profile. 

#3_ParametersDefintion
This file contains the programm used to compute the parameters 

ComputeCorr : Return the correlation length from a segmented profile from SIMULATION or EX-VIVO
ComputePorosity : Return the porosity, the size of pores and the images highlighting the region of interest. The program is used for EX-VIVO.
ComputeInterfaceParameters : Compute the porosity, the Rq and the correlation length of the SIMULATION MAP.	

Old :
Contains programmes and function used to study and test the definition of the parmeters.

#DisplayExVivoParameters
This programme computes the parameter for each of the figure contains in BoneImage and plot a figure with the displayed images and the parameters. 
