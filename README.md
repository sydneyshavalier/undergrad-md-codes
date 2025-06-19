md.f:  The molecular dynamics code I used to perform simulations of molten alkali carbonates. I worked on Lennard-Jones parameterization and worked toward achieving molecular interactions that more closely resembled experimental data.

gr.f: This code calculates the g(r) or radial distribution function between alkali atoms in alkali carbonates.

interface.f: This code creates a box of lithium carbonate molecules to be used in a molecular dynamics simulation.

diffusionk.f: This code calculates the mean squared displacement (MSD) of potassium atoms and plots the resulting value versus time. To find the diffusion coefficient, one only needs to plot the resulting data and multiply the slope of the linear portion of the line by 1/6. While using this code, I would plot the result in xmgrace for further analysis and determination of D.

densityprofile.f: This code plots the densities of both carbon and lithium atoms in lithium carbonate as a function of z. I used histogram binning over a previously recorded molecular trajectory and plotted the results using xmgrace. 

orientation.f: This code tracks the orientations of carbonate molecules in a simulation, and it does so by analyzing the angular orientation of the carbonate plane relative to the Cartesian axes (x, y, z) as a function of z-position in the simulation box. The resulting data is stored in a histogram which displays the orientation versus z when plotted in xmgrace.

dipole.f: This code does calculates the dipole moment of a NaCl system — specifically, it computes the instantaneous total dipole moment vector at each time step from a molecular dynamics trajectory.
