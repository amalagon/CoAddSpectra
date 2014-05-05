This repo will store files that do the Data Analysis for the YMCE experiment.
The pipeline takes in raw FFT spectra, indexed by timestamp and file number,
that are 1D arrays of power density (mW/Hz) vs frequency bin. We fit the spectra empirically, using a rolling average, and subtract the fit to look at the fluctuations of the noise power. The overlapping bins are averaged, and a Lorentzian correction is applied to account for the suppression of axion power off resonance. The final spectra should show excess power in units of the expected axion power vs frequency.
