## Synopsis
<p> This python code was used in the design of a small footprint, high acquisition rate electron spectrometer. The basic design of
an electron spectrometer is to use a region with a known magnetic field to deflect incident electrons onto a detector. Their amount
of deflection will depend on their energy; as a result, each point on the detect corresopnds with a particular electron energy. In
the design code, various parameters including the magnetic field strength, dimensions, and the slit position (where the electrons
enter the detector) were tuned to yield the desired dynamic range. </p>

<p> In the 2D simulation code, electron trajectories are simulated as they enter the magnetic field (mapped out with a hall probe).
Based on their energies and where they impinge on the detector, an energy calibration can be determined. This was calibration was
confirmed experimentally using an electron accelerator. The analysis of the data resulting from this experiment was performed in 
the Absolute_Calibration code. 

<p>
For the detector, a phosphor screen was used in conjunction with a linear CCD. In order to determine the response of this detector
(not known at the high energies used here), the Los Alamos Monte Carlo N-Particle code (MCNP) was used to simulate electron
interactions with a block of material with the same composition as the phosphor screen used. The analysis of the data resulting
from these simulations is performed in the Jupyter notebook titled "Read MCNP Output". In this analysis, it was necessary to know
the effects of photon diffusion within the phosphor material; this was determined in the Photon_Diffusion code. Addtionally, the 
contribution of direct excitation of electron-hole pairs within the CCD from the incident electrons needed to be determined. This
was determined experimentally by blocking regions of the phosphor screen, preventing their fluorescence from reaching the CCD. The
analysis of the resulting data was performed in the Jupyter notebook titled "Lanex_Strip_Test".
</p>

## License
MIT license.
