
# Code used to compute co-firing and co-predcition of units during co-occurring ripples #
original method and results described in Verzhbinsky et al. 2023, PNAS
Some of this is not fully tested, so please reach out with any questions or concerns:
Ilya Verzhbinsky ilya@health.ucsd.edu 

Most of the scripts require three inputs:
- ripple detection run on the LFP data. 
    - Please refer to (https://github.com/iverzh/ripple-detection) for information on how to generate the rippleStats.m output.
- output from unit sorting and classification. This should be saved as a MATLAB cell array variable named 'units' with the following formatting:
    - N x 3 where N is the number of units 
    - Dim 1 indicates the channel the unit was detected on, or where the unit waveform has the largest amplitude.
    - Dim 2 is an array of the spike times in seconds. 
    - Dim 3 is a string that indicates the results of unit classiciation: 
        - 'pyr' for RS putative pyramidal cell
        - 'int' for FS putative interneuron 
        - 'mult' for multiunit 
- a 1 x nSamples logical array named 'StageMask' that indicates the times during the recording that have been staged as a particular behavioral state (e.g. NREM or waking)
