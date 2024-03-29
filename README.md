
# Code used to compute co-firing and co-predcition of units during co-occurring ripples #
Original method and results described in Verzhbinsky et al. 2023, PNAS.
Some of this is not fully tested, so please reach out with any questions or concerns:
Ilya Verzhbinsky ilya@health.ucsd.edu 

Most of the scripts require three inputs:
- ripple detection run on the LFP data. 
    - Please refer to (https://github.com/iverzh/ripple-detection) for information on how to generate the rippleStats.m output.
- output from unit sorting and classification. This should be saved as a MATLAB cell array variable named ``units`` with the following formatting:
    - N x 3 where N is the number of units 
    - Dim 1 indicates the channel the unit was detected on, or where the unit waveform has the largest amplitude.
    - Dim 2 is an array of the spike times in seconds. 
    - Dim 3 is a string that indicates the results of unit classiciation: 
        - 'pyr' for RS putative pyramidal cell
        - 'int' for FS putative interneuron 
        - 'mult' for multiunit 
- a 1 x nSamples logical array named ``StageMask`` that indicates the times during the recording that have been staged as a particular behavioral state (e.g. NREM or waking)
- a nChannel x nSamples matrix called ``RBphaseAll`` that  contains the phase of the ripple bandpass in your LFP data. 


## Computing pairwise unit prediction ##
### ComputeCoRippleCCGs.m
inputs:
- ``units`` cell array
- rippleStats.mat ripple detection output
 - ``StageMask``

output:
 - unit cross-correlograms during three conditions
    - co-rippling periods 
    - periods absent of rippling
    - periods where only one of the cell locations has a detected ripple

### PreprocessPrediction.m
inputs:
- ``units`` cell array
- rippleStats.mat ripple detection output
- the output folder containing all resutls from *ComputeCoRippleCCGs.m*
- RBphaseAll ripple band phase matrix. 
output:
- ``rippStats`` structure containing preprocessed data during co-ripples and ``nullStats`` containing the same for control condition. The most relevant fields of the output structure are:
    - times: contains the relative firing times of unit A and B.
    - trials: the trial ID during which the time in times occurred.
    - APar: the time during the recording that the A unit fired

### ComputePredictability.m
inputs:
- preprocessed ``rippStats`` and ``nullStats`` from *PreprocessPrediction.m*
- ``units`` cell array
- RBphaseAll maxtrix. 

output:
- modified ``rippStats`` and ``nullStats`` structures with the following important fields:
    - ``predict``: mean prediction for cell pair A and B.
    - ``FilterOu``t: the full prediction output distribution of all cell B action potentienals fed into the A-B coupling filter.
    - ``Filter``: an example of the A-B coupling filter constructed from leaving out the last B unit AP.
    - ``Y``: only relevant in ``nullStats``. The trials used to match the number of action potientials between null and co-ripple conditions.
    - ``firingRate``: firing rate of the B unit during the co-ripple or null periods. 
    - ``CellType``: the type of pairwise interaction for A and B (e.g., 'pyr->int')
    - ``PhaseA``: the rippleband phase in location A when unit B fires.
    - ``PhaseB``: the rippleband phase in location B when unit B fires.
    - ``rippPhaseUnitA``: the rippleband phase in location A when unit A fires.

### ComputePredictability.m
All inputs and outputs are the same as *ComputePredictability.m*.
This script computes prediction during a control period that attempts to match the firing rate of unit B for each B->A unit pair, and not the total number of action potentials as is done in *ComputePredictability.m*.

All of the outputs of the predicition code used in Verzhbinsky et al. 2023, PNAS can be found at:
https://zenodo.org/records/10183659

## Figures and analysis ##
These scripts use the outputs deposited in the *out/* directory in the repository https://zenodo.org/records/10183659.

### CoRippleFigures.m
Produces the co-ripple plots in Fig 2.

### CoFireFigures.m
Produces all co firing panels in Fig 3.

### PredictabilityFigures.m
Produces all prediction panels in Fig 4.

### PredictabilityPhaseFigures.m
Produces all prediction phase panels in Fig 5.

### predictStatistics
inputs:
    - ``units`` cell array
    - processed ``rippStats`` and ``nullStats`` from *ComputePredictability.m*

This computes many of the predicion statistics reported in Verzhbinsky et al PNAS 2023

### rippleStatistics
inputs:
    - ``units`` cell array
    - processed ``rippStats`` and ``nullStats`` from *ComputePredictability.m*
    - ``StageMask`` 

This computes many of the ripple-unit (e.g., changes in unit firing rate, ripple-unit phase locking, ... etc.) and ripple characteristics (e.g., ripple densities, phase locking, ... etc) reported in Verzhbinsky et al PNAS 2023

