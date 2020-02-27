# Simple-Fluid-Sim

This is a code sample from my research repo that implements a simple fluid simulation solver. 
The same code does both 2D and 3D simulation: dimensionality is passed in as a template paramter.
- array.h: a N-dimensional array template, equipped with basic arithmetic and iterating tools;
- level_set*.h: implement N-d levelset, with boolean operations, redistancing methods and curvature solve method;
- velocity_field.h: stores N-d velocity field on a staggered grid. With advection methods for scalar fields adn vector fields;
- simulator.h: stores the macroscopic config for the simulation; drives the loop; save/load checkpoints;
- simple.h: implements step/projection methods to be called by Simulator class.

The rest of the repo is constantly modified and mostly experimental. Files above are just a few that has stabilized over time.
