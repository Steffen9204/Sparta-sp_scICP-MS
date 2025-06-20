# Sparta - Single Particle Analysis & Reliable Tracking Algorithm
<b>Suitability</b>: Python post-processing algorithm for single particle / single cell ICP-MS (<b>quadrupole based</b>) suitable for <b>Thermo and Agilent.</b>
<b>Note:</b> The code might be also <b>suitable</b> for <b>other ICP-MS manufacturers</b> if the output is provided as .csv (e.g from Perkin Elmer/Analytik Jena) but the input notebook must be adjusted to the format needed.

<b>Authors Note:</b> We tried to explain how to operate our code as transparent and understandable as possible and hope that our 'ReadMe' helps you to use our Python algorithm. If there are questions please leave comments that we can modify our instructions. Thank you for your support!

## 1)  General Information
* The Python algorithm is created using <b>Anaconda Navigator (2.6.0)</b> and JupiterLab (2.2.6). We suggest to download this software. Otherwise a python script is also available (.py) but is not tested from our site.<br>
* The <b>anaconda environment</b> is provided in the main path <b>('2024_10_base_clone.yaml')</b> and can be integrated in anaconda navigator: Environments --> Import --> Choose File<br>
* After that, it can be chosen as default: File --> Preferences --> Default conda environment --> choose environment. <br>
* Now install the right <b>JupyterLab Version (2.2.6)</b> and 'Launch'.<br>

## 2) Create Folders and insert files
* Create for each sample one folder --> Add <b>ONE csv. file (raw data of the measurement) </b> and the <b>code file </b> (!IMPORTANT: Make sure only one .csv is in the same folder with the code file!). <br>
* <b>Open the code</b> file with e.g. JupiterLab.

## 3) Operate the code (e.g. calculation of transport efficiencies)

### 3.1) Reference material code (example Au)
* Use the reference material code (example Au) to calculate the transport efficiencies (TE) and/or interpret Au particle data. <br>
* Suggested to measure three reference material technical replicates for the transport efficiency and take the average of the TE results. <br>
* <b>In the code</b>: Always when you find <b>Markups with '!!!' ---> User action is required</b> (e.g. second 'notebook': all data of the measurement have to be inserted). <br>
* The <b>intercept and slope (response)</b> is usually taken from a <b>linear ionic calibration</b> (e.g. done with Microsoft Excel or Origin Lab etc.). <br>
* <b>'te'</b> (transport efficiency) here can be for now <b>roughly estimated</b> (as absolute value; NOT as '%') and then <b>updated</b> after the result is there from e.g. three replicates in the code for analyte interpretation (and this code if aiming to investigate Au particles). <br>
* <b>Run the whole code</b>. In case an <b>error</b> occurs, <b>adapt the user interaction fields</b> marked with '!!!'.<br>
* <b>Note:</b> It might be necessary to <b>convert your .csv files to the UTF-8 format</b>, if an error is occuring in the code!
* To decide which <b>method</b> and <b>factor 'k'</b> should be chosen to estimate the transport efficiency, the <b>graphs</b> 'Plot all data and show PDT to choose the best PDT (Gaussian)' & '# Plot all data and show PDT to choose the best PDT (Poisson)' might help. Usually the <b>Gaussian method between µ+4-7SD is taken (!!!to estimate the transport efficiency!!!) in combination with the Particle number method</b> but might vary depending on your particle material used. <br>
* If <b>dwell times < 2 ms</b> are used, a <b>peak integration</b> is applied for both methods to combine multiple detections within one particle event.
* The ionic concentration, which is calculated assuming that all data below PDT are dissolved/ionic signals is calculated. The <b>mode of all data < PDT </b> and the slope & intercept lead to the final <b>ionic concentration (µg/L)</b>, which does NOT include the dilution factor, for each method. <br>
* In the Notebooks <b>'Remove outliers [>µ+kSD from all particles] (high signals e.g. agglomerates) (Gaussian/Poisson)'</b> the <b>factor 'i'</b> might be <b>adapted</b> from default '3' to any higher number to become less sensitive and detect less outliers (in case your particles are detected as outliers, sometimes the case for high background elements (e.g. Si)' This can be monitored in the histogram graphs below. <br>
* The <b>Particle Number Concentration (PNC)</b> shown here is calculated from <b>all events above PDT</b>. <br>
* The below <b>histograms give an overview about background, particles (used for the following calculations) and outliers</b>; Adapt the above factor 'i' if you think that your particle populations are marked as outliers. <br>
* The <b>transport (nebulisation) efficiency</b> is calculated via particle number and particle size method and saved in the final Excel output for all cases. As explained above, take the <b>most suiting method</b> for your materials and calculate an <b>average</b> of e.g. three replicates (always the same method). <br>
* Next, <b>masses, sizes and each limit of detection</b> are calculated and also stored in the Excel output. <br>
* The <b>next interaction</b> might be nessesary in the <b>mass and size histograms (marked with '!!!').</b> <br>
* You can adapt the <b>expected peak number</b> (depending on homogeneous/heterogeneous particle distributions). <b>Note:</b> Sometimes (e.g. if there is leftover background to the left) a higher 'n' suits to a better fit. This has to be tested and can be monitored in the graph(s) below. As described in the Markups: line 19 "(max(h[1]))+"NUMBER"*(max(h[1]))(max(h[1]))+100*(max(h[1]))"; the <b>'Number' can be adapted</b> to adjust the <b>peak form and height</b> (usually between 0.01 and 100). <br>
* In case the <b>peak fitting is chosen</b> as final result (for e.g. heterogeneous materials or obvious background to the left), the <b>events unter the peak can be extracted</b> as well as their <b>particle number concentration</b> (save in the final Excel output file which is explained under <b>'4.'</b>). <br>
* Next the <b>mass concentrations (mg/L)</b> are calculated using the sum of all events and the average times the number of events which is in some case slightly different. We <b>recommend</b> the use of the <b>sum</b>, as this represents all particle events detected. <br>
* Finally, all <b>data of interest</b> are <b>stored and exported</b> alongside the output <b>figures</b> in an <b>Excel file (named with 'output')</b>, which is created in the same folder as you run your code. <br>

### 3.2) All other materials and elements (example SiO<sub>2</sub>)
* Use the <b>analyte code (example SiO<sub>2</sub>)</b> to interpret your <b>particle sample(s)</b> in the folder 'Analytes (all other elements)/Example SiO2'. Choose either Agilent or Thermo (depending on your instrument used). <br>
* <b>All steps explained in '3.1'</b> also apply here with the only difference that the <b>final 'te' calculated</b> should be <b>filled in</b> here (in the second notenbook with all the data). The following steps also apply for the reference materials (e.g. Au).<br>

## 4) Excel output file
* <b>Tab: 'Baseline correction'</b>: In case of a <b>baseline drift</b> the baseline may have been corrected. The <b>ionic mode estimation (background)</b> gives you a <b>first estimation</b> if the <b>Gaussian or Poisson method</b> should be <b>chosen</b> to estimate the <b>Particle Detection Threshold (PDT)</b>. Lockwood et al. 2020 (SPCal) suggested to use the Poisson method for <b>low background elements <= 5 counts otherwise</b> the <b>iterative Gaussian method</b> may be chosen (e.g. for dwell times of 0.1 ms --> 50,000 counts per seconds). <br>
* <b>Tabs: 'PDT_Gaussian/Posisson'</b>: show the <b>PDTs (and critical values (LCs))</b> estimated as a second and <b>most important</b> estimate to choose the <b>most suitable method</b> for final data interpretation and for if the Gaussian method is chosen also factor 'k'. <br>
* <b>Tab: 'TE'</b>: !!! Only appears for the <b>reference material (e.g. Au)</b>!!!. <br>
* <b>Tabs: 'Raw_Gaussian/Poisson'</b>: show the final <b>particle raw data above PDT</b> for each method.
* <b>Tabs: 'Number_conc_Gaussian/Poisson'</b>: represent the <b>particle number concentration [nanoparticle (NP)/mL)] (including dilution factor) above PDT (!!!NO Peak-fitting tool!!!)</b> for each method. <br><br>
* <b>Tabs: 'Masses_Gaussian/Poisson'</b>: represent <b>all particle masses above PDT (!!!NO Peak-fitting tool!!!)</b> for each method. <br>
* <b>Tabs: 'Masses_Gaussian/Poisson_Peaks'</b>: represent <b>extracted particle masses above PDT (!!!Peak-fitting tool!!!)</b> for each method. In the graphs below, the Peaks can be monitored and connected to the peak captions for each method. If n(Peaks) < 4, the column is filled with a '0'.<br>
* <b>Tabs: 'Masses_Summary_Gaussian/Poisson'</b>: show a <b>mean/average, median, SD, Mode of all particle masses above PDT (!!!NO Peak-fitting tool!!!)</b> for each method. Peak 0-3 give you the <b>Gaussian peak-fitting</b> (this has nothing to do with the Gaussian method to estimate the PDT!) <b>average +/- SD (k = 1); the number of particle events</b> considered for this and the <b>particle number concentration under each extracted peak</b> (including the dilution factor). <br>
* <b>Tabs: 'Masses_LOD_Gaussian/Poisson'</b>: represent the <b>limit of detections for particle masses</b> for each method. <br><br>
* <b>Tabs: 'Sizes_Gaussian/Poisson'</b>: represent all <b>particle sizes above PDT (!!!NO Peak-fitting tool!!!)</b> for each method. <br>
* <b>Tabs: 'Sizes_Gaussian/Poisson_Peaks'</b>: represent <b>extracted particle sizes above PDT (!!!Peak-fitting tool!!!)</b> for each method. In the graphs below, the Peaks can be monitored and connected to the peak captions for each method. If n(Peaks) < 4, the column is filled with a '0'.<br>
* <b>Tabs: 'Sizes_Summary_Gaussian/Poisson'</b>: show a <b>mean/average, median, SD, Mode of all particle sizes above PDT (!!!NO Peak-fitting tool!!!)</b> for each method. Peak 0-3 give you the <b>Gaussian peak-fitting</b> (this has nothing to do with the Gaussian method to estimate the PDT!) <b>average +/- SD (k = 1); the number of particle events</b> considered for this and the <b>particle number concentration under each extracted peak</b> (including the dilution factor). <br>
* <b>Tabs: 'Sizes_LOD_Gaussian/Poisson'</b>: represent the <b>limit of detections for particle sizes</b> for each method. <br><br>
* <b>Tabs: 'Mass_conc_Gaussian/Poisson'</b>: represent the <b>mass concentrations (mg/L)</b> calculated as sum for all masses as well as using the masse average multiplicated with the event number for each method. <br>
* <b>Tabs: 'Ionic_conc_Gaussian/Poisson'</b>: represent the <b>ionic concentration</b> which is calculated assuming that all data below PDT are dissolved/ionic signals. The <b>mode of all data < PDT </b> and the slope & intercept lead to the final <b>ionic concentration (µg/L)</b>, which does NOT include the dilution factor, for each method. <br>
