# Farouki space PH curve in Matlab
 This code can be used to compute and draw Pythagorean-hodograph(PH) curves in 3 dimensional space with Farouki's algorithm[^1].

- If the initial/final positions,the initial/final unit tangent vectors, and the arc length of the curve is given, this code computes the corresponding PH curve.

- Additional parameters ($\psi_0,\psi_2$) can be specified if the user wants. The default values are $(\psi_0,\psi_2)=(\pi,0)$.

See the example file (main_example.m) for the usage information.

[^1]:  R. T. Farouki, “Existence of Pythagorean-hodograph quintic interpolants to spatial G 1 Hermite data with prescribed arc lengths,” J. Symb. Comput., vol. 95, pp. 202–216, Nov. 2019.