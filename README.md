# Introduction
This repository is used to store the content discussed in the paper 'Exploration of shock-droplet interaction based on high-fidelity simulation and improved theoretical model
', including the software used, theoretical analysis models, and data analysis Excel.

# Software
The algorithm we use is from MFC. MFC is an open-source tool for solving multi-component, multi-phase, and bubbly compressible flows. It is capable of efficiently solving a wide range of flows, including droplet atomization, shock-bubble interaction, and bubble dynamics (Bryngelson, 2020). Its GitHub code is located at: 'https://github.com/MFlowCode/MFC.git'.

S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2021) Computer Physics Communications 4655, 107396 

@article{Bryngelson_2021,
  title = {{MFC: A}n open-source high-order multi-component, multi-phase, and multi-scale compressible flow solver},
  author = {Spencer H. Bryngelson and Kevin Schmidmayer and Vedran Coralic and Jomela C. Meng and Kazuki Maeda and Tim Colonius},
  journal = {Computer Physics Communications},
  doi = {10.1016/j.cpc.2020.107396},
  year = {2021},
  pages = {107396},
}

# Data
The input file used for calculation is placed in the caltivation directory. 'Mach2.4D22' represents the liquid phase as water, with a shock Mach number of 2.4, and an initial droplet diameter of 22mm. '
'Mach2.4D22_N' represents the liquid phase as n-hexane and the other setup is the same as above.

# Analysis
Theoretical analysis model based on ray tracing method in Python file, the 'pre_ min' file records the pressure changes at various conditions. The reciprocal columns of the data are negative pressure value, negative pressure point x position, and y position, respectively. The 'Mach stem' file records the Mach stem's trajectory obtained from post-processing.

