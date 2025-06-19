md.f:  The molecular dynamics code I used to perform simulations of molten alkali carbonates. I worked on Lennard-Jones parameterization and worked toward achieving molecular interactions that more closely resembled experimental data.

gr.f: This code calculates the g(r) or radial distribution function between alkali atoms in alkali carbonates.

interface.f: This code creates a box of lithium carbonate molecules to be used in a molecular dynamics simulation.

diffusionk.f: This code calculates the mean squared displacement (MSD) of potassium atoms and plots the resulting value versus time. To find the diffusion coefficient, one only needs to plot the resulting data and multiply the slope of the linear portion of the line by 1/6.
