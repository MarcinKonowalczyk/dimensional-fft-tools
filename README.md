## dimensional-fft-tools
This is a MATLAB toolbox for one-dimensional fft operations on multidimensional datasets.

A part of this toolbox - dfun is already published as an independent function [on Mathworks](https://uk.mathworks.com/matlabcentral/fileexchange/63686-dfun). I will maintain it there as a latest stable release. This repository will have the bleeding-edge version in it.

### About
Large datasets (arrays) are, and always will be, an inherent part of physical sciences.
Each dimension represents a different dimension of a dataset while the values themselves are the variable of interest.

For example, a spectroscopist might record the absorbance of the sample of interest (*the value*) as a function of wavelength of the incident light (**1D data**), sample's temperature (**2D**), and time, in minutes, after adding a digestive enzyme to it (**3D**).
Now, lets say that they find themselves being unable to resolve the signal on interest (*the change in the rate of the digestion of the sample across a certain temperature range*) due to noise.
The simplest technical improvement they therefore employ is that for each measurement they record a burst of spectra (**4D**) and turn the light source on and off at a certain frequency (say 100Hz).
This will mean that the signal of interest will oscillate at that frequency and all other special components can be rejected (*forgetting for a moment what actually can happen to a light source is modulated with a large-amplitude square wave*).
The issue they then run into is that their data is a 4D array of absorbance values which they need to perform analysis on - for each 1D slice though the 4th dimension i.e. the burst of spectra, they need to window the data, take the `fft`, select the correct subset thereof, subtract background, and finally integrate the peak of interest.
This adds up to quite a heavy operation which, furthermore, needs to be nested in 4 `for`-loops indexing though all the dimensions. This toolbox aims to help with exactly these kind of tasks.

### ToDo's
 - Go through all the WIP notes and put them into this list
 - Add more toolbox description to this doc
 - Figure out how to avoid using sliceDone matrix. this should improve the performance (<- Major restructuring needed)
 - Add more examples to `dfun` documentation
 - Test `dfun` plotting for different modes and `@plotfuns`
 - Add `@fun` output size > vector to `dfun` for numeric and logical outputs
 - Make sure that if `@fun` output is a string, `dfun` works correctly
 - Allow fof `dfun` & `partcell` to operate on subsets of the dimentions; e.g. 2:5 of length(10) dimention. ?? (Not sure if needed/useful)