# Garcia_et_al_2025

Data analysis codes for "Noncanonical Short-Latency Auditory Pathway Directly Activates Deep Cortical Layers."
Michellee M. Garcia, Amber M Kline, Koun Onodera, Hiroaki Tsukano, Pranathi R. Dandu, Hailey C. Acosta, Michael R. Kasten, Paul B. Manis, Hiroyuki K. Kato.
Nature Communications 2025

Updated 2025 May 21th by Hiroyuki Kato

These are MATLAB scripts to determine linear probe trajectories for analyzing sound response properties of individual brain regions.


## Medial Geniculate Body (MGB) recordings
Note that before running MGB codes, you need to prepare tif files including the brain atlas images.
We used Paxinos and Franklin's The Mouse Brain Atlas 5th edition. We chose the Paxinos Brain Atlas instead of the Allen Brain Atlas for MGN mapping
because MGN boundaries, particularly for MGv, appear smaller in the Allen Brain Atlas compared to those inferred from histological analyses, 
such as calbindin-1 immunostaining. Thus, MGN boundaries defined by the Paxinos Brain Atlas align more consistently with available histological data.
We cannot share Atlas data on GitHub due to copyright restrictions.
Please prepare these files yourself or contact Hiroyuki Kato.

For our MGB analysis, we used Paxinos Mouse Brain Atlas, and included the following range:
AP: Bregma -4.04mm through -2.80mm (11 sections)
ML: Lateral 1.50mm through 2.50mm
DV: Ventral 2.50mm through 4.00mm.

For example, "Bregma-2.80mm_BrainSection.tif" shows a traced brain atlas at Bregma -2.80mm, in the range of lateral 1.5-2.5mm, ventral 2.5-4mm.
Similarly, "Bregma-2.80mm_MGB.tif" shows just the MGB subregions traced from the brain atlas in the same x and y ranges.
If you have access to the digital version of the atlas, using cliping masks to extract that portion of the atlas on Adobe Illustrator is an easy way to create these files. 

Brief descriptions about each script below:

1. create_SectionImage_dataset.m
	This code converts the atlas tif files to SectionImage.mat.

2. extract_MGBoutlines.m
	This code converts the MGB boundaries data to polygons and stores in MGBdrawdata.mat.
	
3. determine_probe_trajectory_MGB.m
	We used this code to determine the stereotaxic coordinates and brain regions of individual linear probe channels in our MGB electrophysiological recordings.
	This code saves layer_info.mat, which stores the information about individual recorded channels. This file will be loaded in all subsequent analyses
	to determine properties of different MGB subregions.
	
	
## Inferior Colliculus (IC) recordings
For our IC analysis, we used AP_histology (https://github.com/petersaj/AP_histology) to align our probe trajectories to the Allen Common Coordinate Framework.
This code generates probe_ccf.mat, which we use to assign each linear probe channel to brain regions.

For visualizing probe trajectories within IC, we used mesh data of individual IC subdivisions extracted from the Allen CCF.
Download structure mask files structure_4.nrrd, structure_811.nrrd, structure_820.nrrd, and structure_828.nrrd, which correspond to IC, ICc, ICd, and ICe, respectively,
from https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_10/ .
We convert these nrrd files to meshdata.

Brief descriptions about each script below:

1. save_surface_mesh.m
	This code converts the nrrd files to surface mesh data and saves "IC_mesh.mat", "ICc_mesh.mat", "ICe_mesh.mat", and "ICd_mesh.mat".
	
2. save_layer_info_IC.m
	We used this code to determine the stereotaxic coordinates and brain regions of individual linear probe channels in our IC electrophysiological recordings.
	This code saves layer_info.mat, which stores the information about individual recorded channels. This file will be loaded in all subsequent analyses
	to determine properties of different IC subregions.

3. draw_probe_trajectories_IC.mapping
	We used this code to display the linear probe trajectories within the IC.
	
	
## Auditory Cortex (AC) recordings
For our AC recordings, we determine layer boundaries using three features.
1. Distribution of spikes is used to determine the brain surface and white matter.
2. Current source density (CSD) analysis during click sound presentations is used to give the first estimation of the top and bottom positions of layer 4.
3. Finally, cross-channel correlation is used to cluster channels into different layers.

Before running the following save_layer_info_auditory_cortex.m, pre-determine the surfaceCh, whitematterCh from the spike data.
CSD signals were derived using CSD.m available on MATLAB Central File Exchange at https://www.mathworks.com/matlabcentral/fileexchange/69399-current-source-density-csd.
Layer 4 (L4) was estimated as the channels in the middle layer with short-latency strong responses.
We used LFP data downsampled to 250Hz.

Brief descriptions about each script below:

1. save_layer_ifo_auditory_cortex
	Cross-channel correlation is used to determine layer boundaries.
	This method typically segregates channels into four groups: layer 1, layer 2-4, layer 5, and layer 6. 
	ref: Senzai et al., Neuron 2019.  https://pubmed.ncbi.nlm.nih.gov/30635232/
	We determine the L3-L4 border based on the L4-top channels in CSD signals, and L4-L5 border as the mean of CSD- and LFP-based boundaries.
	This file will be loaded in all subsequent analyses to determine properties of cortical layers.
