# iolite4_SmNd_UWA_v4.9_debug
Data Reduction Scheme (DRS) script developed in Python for Iolite 4, used for Sm-Nd isotope analysis at the University of Western Australia (UWA). Help needed with visibility of rho and standard-corrected ratios

## Purpose

This DRS script was created to:

- Calculate **standard-corrected isotope ratios** for Sm-Nd (e.g., `143Nd/144Nd`, `147Sm/144Nd`)
- Output a **visible correlation coefficient (rho)** between `log10(StdCorr_147Sm/144Nd)` and `log10(StdCorr_143Nd/144Nd)`
- Improve **spline validation** 
- Clean and streamline outputs for **publication-ready** use

---

## Development & Supervision

- **Code implementation**: Sandra Villacorta  
- **Geochronology supervision**: Chris Fisher (UWA-CET)  
- **Inspiration from**: the original `Sm_Nd_DHF` script by Vitor Barrote, Chris Fisher, and Joe Petrus  

Chris Fisher provided the scientific guidance, geochronology expertise, and reviewed the results, while Sandra implemented the Python logic for Iolite 4.

---

## Scientific References

- Fisher et al., 2020  
- Barrote et al., 2021 – [DOI: 10.5281/zenodo.5512126](https://doi.org/10.5281/zenodo.5512126)

---

## Main Features
- `rho` output is:
  - Registered as an `AssociatedResult` (visible per selection)
  - Plotted as a flat time series (for visibility in output panel)
- Spline-based standard correction for `143Nd/144Nd` and `147Sm/144Nd`

---

## Seeking Input: Result Visibility

Ongoing testing and validation are being conducted to refine the display of output results in Iolite. An observed issue is the inconsistent appearance of standard-corrected ratios (StdCorr_147Sm_144, StdCorr_143Nd_144) and the rho correlation within the Results section, despite their proper AssociatedResult registration. We welcome insights and debugging assistance, especially from the Iolite development team and those with custom DRS framework experience.

---

## Contact

- **Sandra Villacorta** – villacortasp@gmail.com  
- **Chris Fisher** – chris.fisher@uwa.edu.au  

Feel free to open an issue if you have feedback or questions.
