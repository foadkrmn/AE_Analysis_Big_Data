# AE_Analysis_Big_Data
Acoustic Emission signal analysis of Carbon Fiber Reinforced Polymers under fatigue loading

The codes in this respitory are designed to analyze AE signals generated during a fatigue test on CFRP specimen. Two different sensors are used to record data and a typical test generated a total of ~4,000,000 signals. Each signal is described by 16 descriptors and a waveforms with ~15,000 data points. The goal of these codes are:
1 - Synchornize the signal generation time with load time and assign a cycle number to a signal.
2 - Filter out signals that have not recorded in both sensors.
3 - Search and find the waveform of a signal based on time of recording.
4 - Import waveforms and measure the information entropy of the signals.

The main issue in writing and running these codes are managing the large volume of data. This experimental campaign consisted of 10 tests, which resulted in performing all the aformentioned steps for a total of nearly ~40,000,000 signals. 
