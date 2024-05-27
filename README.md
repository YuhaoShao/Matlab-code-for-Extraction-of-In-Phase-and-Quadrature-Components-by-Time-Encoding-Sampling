# Matlab-code-for-Extraction-of-In-Phase-and-Quadrature-Components-by-Time-Encoding-Sampling
The repository contains key MATLAB codes for reproducing the figures from the paper " Extraction of In-Phase and Quadrature Components by Time-Encoding Sampling," authored by Y. H. Shao, S. Y. Chen, H. Z. Yang, F. Xi, H. Hong and Z. Liu. The paper can be accessed online at ***********预印版的网址*************
## 1.	Notes for codes
### Sampling and Reconstruction Codes
BL-TEM.m：for the sampling and reconstruction of the bandlimited signals by the bandlimited-TEM and the bandlimited-POCS method.

BP-TEM.m：for the sampling and reconstruction of the bandpass signals by the bandpass-TEM and the POCS method.

IQ-TEM.m: for the sampling of the bandpass signals using the bandpass TEM and the reconstruction of I/Q components using the alternating POCS method.


### Performance Simulation Codes
SNDRWithoutNoise.m: for the validation of Theorem 1 and Theorem 2.

SNDRBandpassNoise.m: for simulating the reconstruction performance in bandpass noise environments.

SNDRWhiteNoise.m: for simulating the reconstruction performance in white Gaussian noise environments.

SNDRUnderTimingQuanti.m: for simulating the reconstruction performance under timing quantization.

The simulation example and TEM parameters are given in the paper and are specified in the codes. You may replace the signal and TEM parameters for your purposes.

## 2.	Generation of figures in the paper
Figs.4 (a) (curves “by BP-TEM” and “Original”) and (b) are produced by running BP-TEM.
Figs.4 (a) (curves “by BL-TEM”) and (c) are produced by running BL-TEM.

Fig.5 is obtained by running IQ-TEM.m and BL-TEM.m. Running IQ-TEM.m yields the reconstructed I/Q signals (curves “by BP-TEM” in Figs. 5 (a) and (b)). Running the BL-TEM.m yields the reconstructed signals (curves “by BL-TEM” and “Original” in Figs. 5 (a) and (b)) and the outputs of the integrator (Figs.5 (c) and (d)) by replacing the input signals as I/Q components.

Fig.6 is obtained by running SNDRWithoutNoise.m. For a given threshold (cita), running this code outputs a set of SNDR values. Figs.6 (a), (b) and (c) are produced with the threshold (cita) values of 1/120, 1/240, and 1/360, respectively.

Fig.7 (a) is obtained by running SNDRWhiteNoise.m. For a given carrier frequency, a set of SNDR values are outputted. Fig.7 (a) are produced with the carrier frequencies(f0) of 150, 600, 1050 and 1500, respectively.
Fig.7 (b) is obtained by running SNDRBandpassNoise.m. For a given carrier frequency, a set of SNDR values are outputted. Fig.7 (b) are produced with the carrier frequencies (f0) of 150, 600, 1050 and 1500, respectively.

Fig.8 is produced by running SNDRUnderTimingQuanti.m.
