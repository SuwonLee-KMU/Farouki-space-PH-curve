# Farouki space PH curve in Matlab
 This code can be used to compute and draw Pythagorean-hodograph(PH) curves in 3 dimensional space with Farouki's algorithm[^1].

- If the initial/final positions,the initial/final unit tangent vectors, and the arc length of the curve is given, this code computes the corresponding PH curve.
- Additional parameters (psi_0,psi_2) can be specified if the user wants. The default values are (psi_0,psi_2)=(pi,0).

See the example file (main_example.m) for the usage information.



## Class

### spacePH.m

This class represents a PH curve. 
Parameters determining the curve are contained in the properties of this class.
Once the properties are set, transient properties such as `spacePH.controlPoints`  are automatically updated.

#### Optimizing curve constructing parameters

There are infinite number of space PH curves satisfying constraints given by G1 Hermite data(positions and directions), and arc length constraints.
There are two additional design parameters psi0 and psi2 which determines the shape of the unique space PH curve.

Users can obtain the optimal parameters using `spacePH.getOptimalPsi()`.
The optimization cost function is curvature minimization[^2].



### paintPH.m

This class is for visualizing the curve.



### pointPH.m

This class evaluates a PH curve represented as `spacePH` class instance.
When allocated, the evaluating point of curve parameter is 0 by default.
Users can set the curve parameter and evaluation is conducted automatically for transient properties including `pointPH.paremetricSpeed`,`pointPH.curvature`, and `pointPH.torsion`[^3].

For the computations of Bernstein polynomials and its derivatives, [^4] is referenced.



## References

[^1]:  R. T. Farouki, “Existence of Pythagorean-hodograph quintic interpolants to spatial G 1 Hermite data with prescribed arc lengths,” J. Symb. Comput., vol. 95, pp. 202–216, Nov. 2019.
[^2]: R. T. Farouki, C. Giannelli, C. Manni, and A. Sestini, “Identification of spatial PH quintic Hermite interpolants with near-optimal shape measures,” Comput. Aided Geom. Des., vol. 25, no. 4–5, pp. 274–297, May 2008.
[^3]:  R. T. Farouki, C. Giannelli, and A. Sestini, “Helical polynomial curves and double Pythagorean hodographs I. Quaternion and Hopf map representations,” J. Symb. Comput., vol. 44, no. 2, pp. 161–179, Feb. 2009.
[^4]: E. H. Doha, A. H. Bhrawy, and M. A. Saker, “On the Derivatives of Bernstein Polynomials: An Application for the Solution of High Even-Order Differential Equations,” Bound. Value Probl., vol. 2011, pp. 1–16, 2011.