# dangling-orient-clusters
C++ code to analyse the quasi-liquid-layer of ice


Description: Code to analyze lammps binary trajectories. Read in DCD file, cell data, q3 bond order data
 * Code will:
 *       A) build a neighbour list of oxygen atoms and determine the dangling molecules.
 *       B) calculate orientation of molecules from the Euler angles. Ref. frame is xyz.
 *       C) determine cluster sizes of water-like and ice-like molecules used the Depth-first search algorithm.
 *       D) determine the QLL/Ice interface based distribution of coordinates of largest ice-like and water-like clusters.
 *       E) calculate heatmaps of q3, orientation, density and euler angles. Heatmaps calculated for slices along z-axis
 *       F) produce density, q3, orientation profiles.
 *       G) calculate collective statistics such as q3/orientation of dangling mols, ice qll molecules
 *       H) determine number of ice-like and liquid-like molecules
 *
 *   q3 refers to the 3rd order Steinhard parameter and must be calculated beforehand and inputed
