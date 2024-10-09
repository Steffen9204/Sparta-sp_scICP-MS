# Python-algorithm-sp/scICP-MS
Suitability: Python post-processing algorithm for single particle / single cell ICP-MS (quadrupole based) suitable for Thermo and Agilent

<b>1)  General Information:</b><br>
*The Python algorithm is created using Anaconda Navigator (2.6.0) and JupiterLab (2.2.6). We suggest to download this software. Otherwise a python script is also available (.py) but is not tested from our site.<br>
*The anaconda environment is provided in the main path ('2024_10_base_clone.yaml') and can be integrated in anaconda navigator: Environments --> Import --> Choose File<br>
*After that, it can be chosen as default: File --> Preferences --> Default conda environment --> choose environment<br>
*Now install the right JupiterLab Version (2.2.6) and 'Launch'<br>
<br>
<b>2) Create Folders and insert files:</b><br>
*Create for each sample one folder --> Add ONE csv. file output and the code file (!IMPORTANT: Make sure only one .csv is in the same folder with the code file!) <br>
*Open the code file with e.g. JupiterLab
<br>
<b>3) Calculation of transport efficiencies</b><br>
*Use the reference material code (Example Au) to calculate the transport efficiencies (TE) and/or interpret Au particle data <br>
*Suggested to measure three reference material technical replicates for the transport efficiency and take the average of the TE results <br>
*In the Code: Always when you find Markups with '!!!' ---> User action is required (e.g. second notebook: all data of the measurement have to be inserted) <br>
*The intercept and slope (response) is usually taken from a linear ionic calibration (e.g. done with Microsfot Excel or Origin Lab etc.) <br>
*'te' here can be for now estimated and then updated after the result is there from e.g. three replicates <br>
*Run the whole code <br>
*To decide which method and factor 'k' should be chosen to estimate the transport efficiency, the graphs 'Plot all data and show PDT to choose the best PDT (Gaussian)' & '# Plot all data and show PDT to choose the best PDT (Poisson)' might help. Usually the Gaussian method between µ+4-7SD is taken in combination with the Particle number method but might vary depending on your material used. <br>
*In the Notebooks 'Remove outliers [>µ+kSD from all particles] (high signals e.g. agglomerates) (Gaussian/Poisson)' the factor 'i' might be adapted from default '3' to any higher number to become less sensitive and detect less outliers (in case your particles are detected as outliers, sometimes the case for high background elements (e.g. Si)' This can be monitored in the histogram graphs below <br>
*The Particle Number Concentration (PNC) shown here is calculated from all events above PDT.


