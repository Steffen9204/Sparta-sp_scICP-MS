# 2025\_09\_Change\_Log\_Sparta\_1.2

 **1) Transport efficiency - Poisson method calculation fixed.**

* There was a mistake in the calculation, which is now fixed and gives a reasonable transport efficiency.

 **2) Video scan (as '.gif') through raw data to evaluate correct Particle detection threshold (PDT) for both Gaussian and Poisson methods.**

* Instead of a static image as '.jpg', there is now a video visible in the code, which is saved as '.gif' to evaluate the correct PDT. Speed can be adapted, and the video can be stopped at any point.
* In the final Excel output, there is still a static '.png' as '.gif' cannot be imported into Excel files.

**3) Sparta can now be used if no particles are detected (e.g. for Blanks) .**

* If any lists are empty (e.g. no particles above PDT), Sparta reports something like; 'No data available'. This information is also exported into the final Excel Output file.
* Before (Version 1.0 and 1.1), Sparta reported an error if any list was empty and the code stopped. Now the code continues.
* This was done for following Notebooks:

&nbsp;		Ionic data, Gaussian \& Poisson.

&nbsp;		Outlier discrimination, Gaussian \& Poisson.

&nbsp;		Particle number concentration, Gaussian \& Poisson.

&nbsp;		Histograms, Gaussian \& Poisson.

&nbsp;		Transport efficiency, Gaussian \& Poisson.

&nbsp;		Mass calculations, Gaussian and Poisson.

&nbsp;		Mass histograms, Gaussian and Poisson.

&nbsp;		Size calculations, Gaussian and Poisson.

&nbsp;		Size histograms Gaussian and Poisson.

&nbsp;		Particle mass concentration, Gaussian \& Poisson.

&nbsp;		Save mean function, Gaussian \& Poisson.

**4) Output graphs for Masses and Sizes slightly improved.**

* A title was added in case for some PDTs there is no data that it is clear which graph belongs to which method and threshold criteria.



