# SpikefinderCompetition2017
An algorithm for the spikefinder competition : http://spikefinder.codeneuro.org/

The purpose of the competition is to find an algorithm that is fed with calcium imaging time traces of single neurons and is able to output the most likely temporal trace of spiking probabilities.

This algorithm is extremely simple, basically consisting of one simple derivative step, followed by slight smoothing, thresholding and averaging across a small time window. The only adjustable variable is the delay that is used for the derivative.

Using the training data provided by the spikefinder competition, I realized that this optimal delay can be accurately estimated using the measured kurtosis of the calcium time trace, providing a robust and very hands-off algorithm.

This repository includes a Matlab-file that reads in the data and visualized some of its properties (firstGlanceAtDataset.m, unrelated to the spikefinder competition), then the algorithm itself, applied to either the training and the test dataset. I uploaded a small fraction of the original dataset, such that algorithms can be run without any addition after downloading the full content of the repository.

Peter Rupprecht (February 2017)