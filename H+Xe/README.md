# Quantum Scattering & Fano Resonance Analyzer (H+Xe)

A high-performance computational physics pipeline for analyzing angle-dependent Fano resonance profiles in cold atomic collisions. Originally based on the theoretical framework for H+Kr scattering developed by Koike et al., this codebase has been extensively upgraded to model the **H+Xe system** and features a fully automated, dual-run Python visualization engine.

The suite solves the Schrödinger equation using the **Variable Phase Approach (VPA)** to extract macroscopic scattering phase shifts and differential cross-sections (DCS), and synthesizes 3D resonance profiles using the complex-valued $w(\theta)$ parameter.

---

## 📑 Table of Contents

1. [Overview](https://www.google.com/search?q=%23overview)
2. [Codebase Architecture](https://www.google.com/search?q=%23codebase-architecture)
3. [Key Features](https://www.google.com/search?q=%23key-features)
4. [Modifications from the Original Framework](https://www.google.com/search?q=%23modifications-from-the-original-framework)
5. [Prerequisites & Installation](https://www.google.com/search?q=%23prerequisites--installation)
6. [Usage Guide](https://www.google.com/search?q=%23usage-guide)
7. [Expected Outputs](https://www.google.com/search?q=%23expected-outputs)
8. [References](https://www.google.com/search?q=%23references)

---

## 🔬 Overview

In cold atomic collisions (where the de Broglie wavelength is comparable to the interaction range), atoms can temporarily trap each other behind centrifugal barriers, creating **orbiting shape resonances**. Because different partial waves ($l$) interfere, the resulting resonance profiles (Beutler-Fano profiles) change their shape drastically depending on the observation angle $\theta$.

This codebase computationally models this phenomenon for the **H+Xe** interaction, parameterized by the Toennies et al. (1979) hybrid Lennard-Jones/dispersion potential. It calculates the complex parameter $w(\theta)$ and the Fano $q(\theta)$ parameter to map the topological evolution of the resonance across the complex plane.

---

## 🏗 Codebase Architecture

The project is structured into a classic two-tier scientific computing pipeline:

* **`/vpa/` (Tier 1: The Physics Engine):** Contains the Fortran 90/77 source code.
* `Main.f90`: The driver routine that parses control words (`DCS`, `WPRM`, `ESCAN`).
* `Potential.f90`: Defines the physical interaction (L-J well depth, equilibrium distance, $C_6, C_8, C_{10}$ dispersion coefficients).
* `VPALib.f90`: Core quantum mathematics synthesizing the Clebsch-Gordan interference.
* Libraries: Integrators (`DLSODALIB.f`), Riccati-Bessel functions (`RICBES.f`), and angular momentum algebra (`cgcoef.f`).


* **`/plots/` (Tier 2: The Visualization Pipeline):** Contains the Python master script `generate_final_plots.py`. It automates the Fortran engine, handles complex data parsing, and orchestrates `Matplotlib` and `Gnuplot` to generate publication-quality PDFs and 3D rotating GIFs.

---

## ✨ Key Features

* **Variable Phase Approach (VPA):** Directly integrates the phase shift differential equations without needing to match wavefunctions at the origin.
* **Dual-Run Automation:** The Python script automatically runs the Fortran engine twice—once to calculate overarching macroscopic physics (Phase Shifts, Standard DCS) and once to zoom in precisely on the microscopic Fano resonance.
* **Complex $w(\theta)$ Trajectory Tracking:** Maps the evolution of the Fano resonance as a continuous spiral in the 3D complex plane.
* **Auto-Scaling 3D Visuals:** Dynamically calculates bounding boxes for complex plane data to ensure static PDFs are large and detailed, while 3D animated GIFs remain perfectly stable on a cubic base.
* **Fano Archetype Slicing:** Instead of plotting arbitrary angles, the codebase mathematically scans the $q$-parameter array to plot the exact angles where the resonance acts as a pure peak ($q \to \infty$), pure dip ($q \to 0$), or maximum asymmetry ($q \approx \pm1.5$).

---

## 🚀 Modifications for System Migration: H+Kr to H+Xe

This codebase was explicitly migrated from the original H+Kr framework (Koike et al.) to model the heavier, more polarizable H+Xe system. The specific physical and computational adjustments required for this transition include:

1. **Interaction Potential Overhaul:**
* Replaced the H+Kr Lennard-Jones/dispersion matrix with the Toennies et al. (1979) empirical parameters for H+Xe in `Potential.f90`.
* The H+Xe van der Waals well is significantly deeper ($-7.08 \text{ meV}$ compared to $\approx -5.90 \text{ meV}$ for H+Kr), requiring an updated equilibrium distance ($7.219 \text{ a}_0$) and much larger $C_6, C_8,$ and $C_{10}$ long-range dispersion coefficients.


2. **Targeting High Angular Momentum Resonances:**
* Because of the deeper potential well, the H+Xe system is capable of trapping quasi-bound states against much higher centrifugal barriers than H+Kr.
* We shifted the computational search grid and VPA solver parameters from low angular momentum states (e.g., $l=4$ orbiting resonance at $4.13 \text{ cm}^{-1}$ in H+Kr) to highly excited shape resonances. The modified engine now precisely targets $l=5, 6, 7,$ and $8$, scaling up to collision energies of $E_r \approx 40 \text{ cm}^{-1}$.


3. **Topological Adaptation for $w(\theta)$ and $q(\theta)$:**
* Migrating to higher $l$-states fundamentally altered the mathematical geometry of the Fano interference. The Legendre polynomial $P_l(\cos\theta)$ for the resonant wave now contains up to 8 angular roots (for $l=8$) instead of 4.
* Consequently, the code now handles a much higher topological complexity: the complex $w(\theta)$ parameter trajectories execute up to twice as many infinite excursions (loops) in the complex plane, and the Beutler-Fano $q(\theta)$ profile parameter undergoes rapid, high-frequency asymmetric phase flips across the $0^\circ \to 180^\circ$ scattering range compared to the simpler H+Kr topology.

---

## 💻 Prerequisites & Installation

**Required Dependencies:**

* `gfortran` (GNU Fortran compiler)
* `make`
* `python 3.8+` (with `numpy` and `matplotlib`)
* `gnuplot` (Must be installed and accessible in the system PATH)

**Installation:**

```bash
git clone https://github.com/YourUsername/H-Xe-Fano-Resonance.git
cd H-Xe-Fano-Resonance
pip install numpy matplotlib

```

---

## ⚙️ Usage Guide

You do not need to manually compile the Fortran code. The Python master script handles the build process, execution, and plotting automatically.

1. **Prepare an Input File:** Create or modify a text file (e.g., `VPA_5input.TXT`) in the `/vpa/` directory. It must start with the `WPRM` control word and contain the appropriate background phase shifts for the target resonance.
2. **Run the Master Pipeline:**
```bash
cd plots
python generate_final_plots.py ../vpa/VPA_5input.TXT

```



The script will automatically trigger the Fortran dual-run, parse the outputs, and generate two organized directories:

* `/Global_System_Plots/`: Contains the macroscopic phase shifts, potential curves, and standard DCS surfaces.
* `/l_5_plots/` (or respective $l$): Contains the high-resolution Fano energy profiles, $q$-parameter graphs, complex $w$ trajectories, and rotating GIFs.

---

## 📊 Expected Outputs

Upon successful execution, the target directories will be populated with 14 publication-ready figures, including:

* **`energy_vs_phase_shift.png`**: Standard phase shift tracking showing the $\pi$-radian jumps at resonances.
* **`energy_angle_vs_asydcs.pdf` & `.gif**`: The synthesized 3D interference surface of the Fano resonance.
* **`angle_vs_q_parameter.pdf`**: Fano $q$-parameter mapped against scattering angle showing asymptotes at the Legendre roots.
* **`real_w_imag_w_vs_angle_rotating.gif`**: A smooth, stabilized 3D animation of the complex $w(\theta)$ vector spiraling through the complex plane.
* **`energy_vs_dcs_profiles.pdf`**: The classic Beutler-Fano 2D slices showing the pure peak, dip, and asymmetric profiles mathematically extracted from the data.

---

## 📚 References

1. **Koike, F., et al.** (2026). *Scattering angle dependence of Fano resonance profiles in cold atomic collisions analyzed with the complex valued w parameter.* Indian J Phys.
2. **Toennies, J. P., Welz, W., & Wolf, G.** (1979). *Molecular beam scattering studies of orbiting resonances and the determination of van der Waals potentials for H-Ne, Ar, Kr, and Xe.* The Journal of Chemical Physics, 71(2), 614-642.
3. **Fano, U.** (1961). *Effects of Configuration Interaction on Intensities and Phase Shifts.* Physical Review, 124(6), 1866.

---

*Developed for research in theoretical atomic scattering and cold collision quantum mechanics.*
