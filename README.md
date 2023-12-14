# Introduction
This warehouse is used to store the content discussed in the paper 'Exploration of shock-droplet interaction based on high-fidelity simulation and improved theoretical model
', including the software used, theoretical analysis models, and data analysis Excel.

# Softwave
The algorithm we use is from MFC. MFC is an open-source tool for solving multi-component, multi-phase, and bubbly compressible flows. It is capable of efficiently solving a wide range of flows, including droplet atomization, shock-bubble interaction, and bubble dynamics (Bryngelson, 2020). Its GitHub code is located at: 'https://github.com/MFlowCode/MFC.git'.

# Data
The input file used for calculation is placed in the caltivation directory. 'Mach2.4D22' represents the liquid phase as water, with a shock Mach number of 2.4, and an initial droplet diameter of 22mm. '
'Mach2.4D22_N' represents the liquid phase as n-hexane and the other setup is the same as above.

# Analysis
Theoretical analysis model based on ray tracing method in Python file, the 'pre_ min' file records the pressure changes at various conditions. The reciprocal columns of the data are negative pressure value, negative pressure point x position, and y position, respectively.
