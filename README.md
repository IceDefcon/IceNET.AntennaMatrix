# IceNet.AntennaImpedanceMatrix



The Source code allow to generate Coupling Impedance Matrix for the 20 triangles copper bar Antenna !!!

Valid Matrix is used to approximate complex IQ samples in Electromagnetic Computation Process

Computation use Method of Moments to Estimate the values of Mixed Potential:

  1. Magnetic Vector Potential
  2. Electric Scalar Potential

Mixed Potentials are used to Approximate values for Electric currents induced on the surface of Antenna

  1. The Excitation Matrix (E-Field) is acted on the Coupling Impedance Matrix
  2. Induced Currents can be calculated when E-Field hits Antenna surface
  3. Valid IQ samples can be derived from the Over All Computation !!!

# ---===[ SOURCE CODE ]===---

Compile:

	./build

Execute:

	./Impedance_Matrix

Clean:

	./clean

Execution ---> gives 19x19 Coupling Impedance Matrix 

Matrix is designed for a Dipole Antenna constructed from 20 triangular patches

Coupling Impedance is calculated between the pair of the Test and Source triangular patches

# ---===[ ANTENNA PHYSICS ]===---

Dipole antena is constructed from 20 triangular patches, 10 squares (each contain two equal triangular patches)

Each square is 0.05m x 0.05 such that the antenna will have the dimentions ---> 50x5 cm 

Electromagnetic Wave Frequency = 300Mhz

This gives wavelength of: 1m ---> so the half wave = 0.5m ( which is Antenna dimentions: Half wave dipole )  

Next step is to:

 	1. Compute Exitation Matrix 
 	2. Act Exitation Matrix at Coupling Impedance Matrix 
 	3. Generate Radiation Pattern (Donut for the half wave dipole)


# ICE
