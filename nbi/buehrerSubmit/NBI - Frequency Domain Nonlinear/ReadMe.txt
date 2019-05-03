The files in this directory support frequency domain cancellation of narrowband interference using an adaptive threshold.  The following functions are included:

FreqDomainCancel.m - This function performs the frequency domain cancellation.  The input signal is converted to the frequency domain and compared to a provided threshold.  All values above the threshold are clipped but while keeping the phase consistent.

calculate_threshold.m - This function determines the adaptive threshold for interference mitigation.  Specifically, it divides the frequency range into 10 sub bands, with the assumption that NBI is isolated to one sub band or possibly two if it straddles to sub bands.  The 8 bands with the lowest energy are used to calculate the threshold since they should represent the signal only.

TestNbiFreqDomain.m - This is a script which tests the frequency domain cancellation technique and demonstrates how to call the functions.