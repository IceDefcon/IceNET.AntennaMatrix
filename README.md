# IceNet.AntennaImpedanceMatrix

This is 1.0.0 version of the source code

Modified Source code allow to generate Coupling Impedance Matrix of the custom shape Antenna !!!

Valid Matrix is used to approximate complex IQ samples in Electromagnetic Computation Process

Computation use Method of Moments to Estamate Mixed Potential values:

  1. Magnetic Vector Potential
  2. Electric Scalar Potential

Mixed Potentials are used to Approximate values for Electric currents induced on the surface of Antenna

  1. The Excitation Matrix (E-Field) is acted on the Coupling Impedance Matrix
  2. Induced Currents can be calculated when E-Field hits Antenna surface
  3. Valid IQ samples can be derived from the Over All Computation !!!

# ---===[ SOURCE CODE ]===---

This code allow you to calculate single element of Coupling Impedane Matrix

To compile:

	gfortran -c definitions.f90 functions.f90
	gfortran  definitions.o functions.o main.f90 -o single_impedace_element

To execute:

	./single_impedace_element

Execution should give one complex value from Entire Matrix for a choosen Antenna Surface Shape !!!

	0.3182E-02 +j -0.2432E-02

# ICE
