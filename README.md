# dangling-orient-clusters
## C++ code to analyse the quasi-liquid-layer of ice


**Description:** Code to analyze lammps binary trajectories. 

**Input**: DCD file, cell data, q3 bond order data.
* **Code will:**
  - A) build a neighbour list of oxygen atoms and determine the dangling molecules.
  - B) calculate orientation of molecules from the Euler angles. Ref. frame is xyz.
  - C) determine cluster sizes of water-like and ice-like molecules used the Depth-first search algorithm.
  - D) determine the QLL/Ice interface based distribution of coordinates of largest ice-like and water-like clusters.
  - E) calculate heatmaps of q3, orientation, density and euler angles. Heatmaps calculated for slices along z-axis
  - F) produce density, q3, orientation profiles.
  - G) calculate collective statistics such as q3/orientation of dangling mols, ice qll molecules
  - H) determine number of ice-like and liquid-like molecules
 
 q3 refers to the 3rd order Steinhard parameter and must be calculated beforehand and inputed

### Output for every frame:
 - **cluster_time.txt:** size of the largest water cluster on the top and bottom surface and the size of the largest ice-like clusters within the top and bottom QLL surfaces. 
 - **dang_time.txt:** number of dangling molecules on the top and bottom surfaces and the average dipole angle and q3 value of the dangling molecules. 
 - **interface_time.txt:** position of the QLL/ice interface and estimate of the QLL width. 
 - **noxygens_time.txt:**  Total number of ice and crystal molecules, number of QLL molecules on each surface, and the average dipole angle and q3 of QLL molecules. 
 
### Output averaged over entire trajectory:
 - **ProfileDensity.txt:** density profile along surface normal
 - **ProfileOrient.txt:** dipole angle profile along surface normal
 - **ProfileQ3.txt:** q3 profile along surface normal
 - **ProfileSlice_dipoleTheta.txt:** dipole angle profile for 50 slices along the surface normal
 - **Heatmap_Density.txt:** Density heatmap for 50 slices along the surface normal
 - **Heatmap_DipoleTheta.txt:** Dipole angle heatmap for 50 slices along the surface normal
 - **Heatmap_OrientationThetaGamma.txt:** Theta (dipole) and Gamma heatmaps for 50 slices along the surface normal
 - **Heatmap_Q3.txt:** q3 heatmaps for 50 slices along the surface normal 
