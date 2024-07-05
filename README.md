This code computes the phase function and the isochrons of the Van der Pol model, using forward integration.
The code can be easily adapted to other models with a limit cycle. The method also works with models of higher dimensions.
Note that the use of the "contour" matlab function is delicate with the phase function, due to the singularity between -pi and pi. We rather use the "contourcs" function (Copyright (c)2010, Takeshi Ikuma).

The method is based on the results presented in "A. Mauroy and I. Mezic, On the use of Fourier averages to compute the global isochrons of (quasi)periodic dynamics,
Chaos, vol. 22, no. 3, p. 033112, 2012"

For more information, please email me at "alexandre.mauroy 'at' unamur.be"
