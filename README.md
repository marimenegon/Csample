# C program sample

Code in C to analyse the output of molecular dynamics simulations.

Simulation of bead-spring chains at various pressures. 
At each pressure, it calculates the order parameters that quantitfy 
(1) the orientation of chains, 
(2) the organisation in layers, and 
(3) the hexagonal ordering of nearest neigbours.

Output is presented in three data files: 
in file 1a, there are the time evolution of global quantities;
in file 1b, there are the quantities relative to each particle, for each configuration;
in file 1c, the final configuration is recorded and used as initial configuration for a new simulation.
