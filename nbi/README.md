Narrowband Interference Mitigation  
Author: Jianyuan (Jet) Yu  
Contact: jianyuan@vt.edu  


# Files
* `nbiMain.m` - main file
* `fftThre.m`, `kayEst`, `trainNF`


# Updates
* (20th Feb) fftThre, Kay, TrainNF implemented 
* (17th Feb)

# Overview
* Time domain narrowband interference mitigation (signal of interest is single-carrier signal with linear modulation) - needs data sheet also;  
*  Frequency-domain narrowband interference mitigation (same signal of interest) - needs data sheet also

## sigal of interest
* bandwidth signal, single carrier, OFDM
* RC pulse shaping

## signal of non-interest
* narrowband, here single tone

## method
### 1) time domain
* FIR filter, truncated impulse response 
    * ZF
    * MMSE
    * LMS
    * RLS

### 2) frequency domain
* threshold filter

### 3) narrowband parameter estimation
* [Kay]()

# Further Work
* DS-SS
* UWB
* OFDM



# Related Files
[data sheet (Google Doc)](https://docs.google.com/document/d/1dULYCHmw7NCb9Pp9SOaIqyIlrDgNPNYFvnaRa6fGYws/edit?usp=sharing)