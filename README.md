# Python-algorithm-sp/scICP-MS
<b>Suitability</b>: Python post-processing algorithm for single particle / single cell ICP-MS (<b>quadrupole based</b>) suitable for <b>Thermo and Agilent.</b>

## 1)  General Information
*The Python algorithm is created using <b>Anaconda Navigator (2.6.0)</b> and JupiterLab (2.2.6). We suggest to download this software. Otherwise a python script is also available (.py) but is not tested from our site.<br>
*The <b>anaconda environment</b> is provided in the main path <b>('2024_10_base_clone.yaml')</b> and can be integrated in anaconda navigator: Environments --> Import --> Choose File<br>
*After that, it can be chosen as default: File --> Preferences --> Default conda environment --> choose environment. <br>
*Now install the right <b>JupiterLab Version (2.2.6)</b> and 'Launch'.<br>

## 2) Create Folders and insert files
*Create for each sample one folder --> Add <b>ONE csv. file (raw data of the measurement) </b> and the <b>code file </b> (!IMPORTANT: Make sure only one .csv is in the same folder with the code file!). <br>
*<b>Open the code</b> file with e.g. JupiterLab.

## 3) How to operate the code (e.g. calculation of transport efficiencies)
### 3.1) Reference material code (example Au)
*Use the reference material code (example Au) to calculate the transport efficiencies (TE) and/or interpret Au particle data. <br>
*Suggested to measure three reference material technical replicates for the transport efficiency and take the average of the TE results. <br>
*<b>In the code</b>: Always when you find <b>Markups with '!!!' ---> User action is required</b> (e.g. second 'notebook': all data of the measurement have to be inserted). <br>
*The <b>intercept and slope (response)</b> is usually taken from a <b>linear ionic calibration</b> (e.g. done with Microsoft Excel or Origin Lab etc.). <br>
*<b>'te'</b> (transport efficiency) here can be for now <b>roughly estimated</b> and then <b>updated</b> after the result is there from e.g. three replicates in the code for analyte interpretation (and this code if aiming to investigate Au particles). <br>
*<b>Run the whole code</b>. In case an <b>error</b> occurs, <b>adapt the user interaction fields</b> marked with '!!!'.<br>
*To decide which <b>method</b> and <b>factor 'k'</b> should be chosen to estimate the transport efficiency, the <b>graphs</b> 'Plot all data and show PDT to choose the best PDT (Gaussian)' & '# Plot all data and show PDT to choose the best PDT (Poisson)' might help. Usually the <b>Gaussian method between µ+4-7SD is taken (!!!to estimate the transport efficiency!!!) in combination with the Particle number method</b> but might vary depending on your particle material used. <br>
*In the Notebooks <b>'Remove outliers [>µ+kSD from all particles] (high signals e.g. agglomerates) (Gaussian/Poisson)'</b> the <b>factor 'i'</b> might be <b>adapted</b> from default '3' to any higher number to become less sensitive and detect less outliers (in case your particles are detected as outliers, sometimes the case for high background elements (e.g. Si)' This can be monitored in the histogram graphs below. <br>
*The <b>Particle Number Concentration (PNC)</b> shown here is calculated from <b>all events above PDT</b>. <br>
*The below <b>histograms give an overview about background, particles (used for the following calculations) and outliers</b>; Adapt the above factor 'i' if you think that your particle populations are marked as outliers. <br>
*The <b>transport (nebulisation) efficiency is calculated via particle number and particle size method and saved in the final Excel output for all cases. As explained above, take the most suiting method for your materials and calculate an average of e.g. three replicates (always the same method). <br>

### 3.2) All other materials and elements (example SiO<sub>2</sub>)
*All steps explained in '3.1' also apply here with the only difference the the final 'te' calculated should be chosen here (in the second notenbook with all the data). <br>
*Next, masses, sizes and each limit of detection are calculated and also stored in the Excel output. <br>
*The next interaction might be nessesary in the mass and size histograms (marked with '!!!')
*You can adapt the expected peak number (depending on homogeneous/heterogeneous particle distributions). Note: Sometimes (e.g. if there is leftover background to the left) a higher 'n' suits to a better fit. This has to be tested and can be monitored in the graph(s) below. As described in the Markups: line 19 "(max(h[1]))+"NUMBER"*(max(h[1]))(max(h[1]))+100*(max(h[1]))"; the 'Number' can be adapted to adjust the peak form and height (usually between 0.01 and 100). <br>
*In case the peak fitting is chosen for e.g. heterogeneous materials or obvious background to the left, the events unter the peak can be extracted as well as their particle number concentration (save in the final Excel output file which is explained under '4.'). <br>
*Next the mass concentrations (mg/L) are calculated using the sum of all events and the average times the number of events which is in some case slightly different. We recommend the use of the sum, as this represents all particle events detected. <br>
* Finally, all data of interest are stored and exported alongside the output figures in an Excel file, which is created in the same folder as you run your code. <br>


<br>
Authors Note: The Readme will be completed as soon as possible!


