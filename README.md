# Neutron Flux Visualization 

## About
Understanding the distribution of neutrons inside a reactor core is essential for nuclear engineering tasks. The interaction of these neutrons determine optimal fuel rod placement and geometry. Using [MCNP 6.2](https://rsicc.ornl.gov/codes/ccc/ccc8/ccc-850.html), we can generate a Monte Carlo simulation of the position of neutrons given various fuel rod designs. 

In this example project, we take the point cloud data from MCNP6.2 and generate a scalar field, using [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) to create a geometry. From here, we can use OpenGL to visualize this data. Our visualization technique uses RGB color mapping and height mapping to demonstrate neutron density at various level sets.

## Run Instructions
Requires Visual Studio. To run, load learnply.sln solution file and execute learnply.cpp.

## Acknowledgements
Project was completed as a part of CS-553 Scientific Visualization @ OSU. Template code provided by Dr. Eugene Zhang.
Thanks to Massimo Larsen for providing simulation data and subject matter advice.