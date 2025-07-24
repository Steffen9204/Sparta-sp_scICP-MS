# 2025_07_Change_Log
<b> 1) 'Clear_Cut' function now disabled by default </b>
* Clear_Cut Function to overwrite the iterative Gaussian Particle Detection Throshold, can be enabled at the beginning if using a monodisperse particle sample and the
  (ionic) background is clearly separated from the particle distribution. In this case, it removes potential background artefacts potentiall detected by the statistical    Gaussian method. Often useful for low-background elements such a Au.

<b> 2) 'Outlier removal' function now disabled by default </b>
* For well-characterized particles such as certified reference materials, which have monodisperse distributions, the ‘Outlier removal’ function can be enabled to
  exclude large agglomerates or other artefacts. We suggest to leave this function disabled for unknown or natural samples, as detected agglomerates may in fact 
  represent large particles.
* Further information can be found in the README file.

<b> 3) 'Poisson particle detection threshold (PDT)' was corrected </b>
* There was a gap between the PDT and the first particle in version 1.0. This was corrected subtracting the average of the background (same as for the particle events) 
  from the PDT instead of the LC. In version 1.1, the correct Poisson size and mass PDTs are now shown and exported into the Output Excel file.

 <b> 4) 'Alpha/Beta errors added for the Gaussian peak-fitting/detection method </b>
* Sparta uses a flexible multi-modal Gaussian peak-fitting algorithm that can model up to four peaks, representing either particle distributions or background signals 
  (details, see README file). The algorithm includes both alpha (α) and beta (β) error calculations:
  * α error defines the confidence interval of the peak height and width (set to '0.05' by default).
  * β error reflects the probability of missing a true peak (Type II error), evaluated via a statistical t-test, returning the test's power (1–β).
* Therefore, by default, the fitting uses a 95% confidence level (modifiable via the 'alpha_error' parameter). Peaks are statistically evaluated to determine which are 
  significant and therefore considered as 'real' peaks.
* In the graphs and Excel output, 'bold' marked peaks are within the 95% confidence interval, all other 'not bold' marked peaks are outside this statistical interval.
* Excel output: Tab(s): 'Masses/Sizes_Summary_Gaussian/Poisson' represents e.g. peak 0 of 3SD: Average[+/- α error] +/- SD [+/- α error] fg/nm with power: 1-β for 
  X events and Y particles/mL.

<b> 5) Instrument-specific code files merged (Agilent vs. Thermo) </b>
* Only one code file independent from the instrument used (version 1.0 had instrument-specific code files).
* In the second notebook (beginning of the code), the instrument can be sepcified for the correct output: "instrument = 'Agilent' # Either 'Agilent' or 'Thermo'".
* It is possible to adapt this part in the code for other output forms such as from Perkin Elmer or similar instruments.

<b> 6) Poisson transport efficiency (TE) now exported to the 'Output Excel file' </b>
* If using the TE code file: TE Output can be found in the Excel Tab ' TE_Poisson'.

<b> 7) Limit of detetcion (LOD) [MassHunter] was corrected to Background equivalent mass/diameter [MassHunter] in the code files' </b>
* Background Equivalent Mass [MassHunter] [fg/event].
* Background Equivalent Diameter [MassHunter]  [nm].
