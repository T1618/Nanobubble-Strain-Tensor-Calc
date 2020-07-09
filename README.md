# Nanobubble-Strain-Tensor-Calc
Matlab scripts for calculate the strain from nanobubbles in atomically thin van der Waals materials

05/25/2020
This repository hosts Matlab code for calculating the strain distributions of nanobubbles in atomically thin materials, e.g. graphene and the transition metal dichalcogenides, based on topography measurements made with atomic force microscopes. 
The repository has four files: the first and principal is the function nanobubble_straintensor_solve.m which solves for the strain tensor and Airy stress function based on a given topographic profile and Poisson ratio of bubbled membrane. The second, ExampleStrainCalc.m is a script that shows how to call nanobubble_straintensor_solve.m to calculate the strain tensor for data Fig3AFM.txt. The final is trace.png which shows the sum of the diagonal strain tensor components for the example data using the given parameters in ExampleStrainCalc.m.
More information on this calculation method can be found in the paper: “Facile and quantitative estimation of strain in nanobubbles of arbitrary symmetry in 2D semiconductors verified using hyperspectral nano-optical imaging” currently submitted to the Journal of Chemical Physics. 

For troubleshooting and other questions please contact Tom Darlington at td2583@columbia.edu.

07/08/2020
Paper is now online.  https://doi.org/10.1063/5.0012817
