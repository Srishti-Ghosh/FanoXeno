# H+Xe Cold Collision Scattering Calculation

This directory contains a modified Fortran VPA (Variable Phase Approach) codebase configured for calculating scattering observables for the **H + Xe** system using the hybrid potential model from **Toennies et al. (1979)**.

## Overview

The code solves the radial Schrödinger equation for atom-atom collisions and computes:
- **Phase shifts** vs. collision energy and angular momentum quantum number
- **Angular differential cross sections (DCS)** vs. scattering angle and energy
- **Interaction potential** profile

### Reference
Toennies, J. P., et al. (1979): "The scattering of H, D, He by N₂, O₂, and Xe molecules and their repulsive interactions" *J. Chem. Phys.* **71**, 614.

---

## Quick Start (macOS)

### 1. Install Fortran Compiler

```bash
brew install gcc
```

Verify installation:
```bash
gfortran --version
```

### 2. Build the Executable

```bash
cd /path/to/CodesForColdCollisionProject/vpa
rm -f *.o difelsccs
/opt/homebrew/bin/gfortran \
  -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk \
  -mmacosx-version-min=10.15 \
  -O2 -o difelsccs Main.f90 VPA-Lib.f90 RICBES.f Potential.f90 DLSODALIB.f WKB.f90 cgcoef.f legendre.f90
```

Check the binary:
```bash
file difelsccs
# Should output: Mach-O 64-bit executable arm64
```

### 3. Run the Scattering Calculation

```bash
cd /path/to/CodesForColdCollisionProject/vpa

./difelsccs << 'EOF'
DCS
0,10,
1.0,50.0,200
0.0,180.0,180
EOF
```

**Input parameters:**
- `DCS`: Control mode (Angle Differential Cross Sections)
- `0,10`: Angular momentum range (l_min=0 to l_max=10)
- `1.0,50.0,200`: Energy range (1 to 50 cm⁻¹) with 200 points
- `0.0,180.0,180`: Scattering angle range (0° to 180°) with 180 points

**Output files generated:**
- `potential.txt`: Interaction potential V(r) in atomic units
- `phshift.txt`: Scattering phase shifts η_l(E)
- `angdcs.txt`: Angular differential cross sections dσ/dΩ(E, θ)

### 4. Generate Plots

```bash
cd /path/to/CodesForColdCollisionProject/analysis
python3 plot_h_xe_results.py
```

This creates four publication-quality plots in `plots/`:
- `potential.png`: Interaction potential curve
- `phshift.png`: Phase shifts vs energy for all l values
- `dcs_vs_angle.png`: DCS at selected collision energies
- `dcs_heatmap.png`: 2D contour of DCS(Energy, Angle)

---

## Code Modifications for H+Xe

The following changes were made to the original H+Kr codebase to use H+Xe parameters:

### Modified Files

#### 1. **vpa/Potential.f90** — Potential Parameters
- **Masses** (in electron masses):
  - `amass_H`: 1837.15 (from 1837.36 for Kr)
  - `amass_Xe`: 239332.50 (replaces 152754.0 for Kr)

- **Lennard-Jones Parameters** (r < r_s):
  - `RE1`: 7.219 a₀ (equilibrium distance, was 6.74632)
  - `EPS`: 0.000260 Hartree (well depth = 7.08 meV, was 0.000216821)
  - `rs`: 8.957 a₀ (switching radius, was 9.1481)

- **Dispersion Coefficients** (r ≥ r_s):
  - `c6`: 42.298 Hartree·a₀⁶ (was 28.62)
  - `c8`: 918.37 Hartree·a₀⁸ (was 573.737)
  - `c10`: 27554.64 Hartree·a₀¹⁰ (was 14837.1)

- **Reduced Mass Calculation**:
  - Changed from `RM = amass_H * amass_Kr / (amass_H + amass_Kr)`
  - To: `RM = amass_H * amass_Xe / (amass_H + amass_Xe)`

#### 2. **vpa/Main.f90** — Mass Update
- Updated `DATA amass_Kr` to 239332.50 electron masses (Xe mass)

#### 3. **vpa/arc/wparam.f90** — Mass Updates
- Updated `DATA amass_H` to 1837.15 electron masses
- Updated `DATA amass_Kr` to 239332.50 electron masses (repurposed variable for Xe)

#### 4. **vpa/arc/VPA-main.f90** — Mass Update
- Updated `DATA amass_Kr` to 239332.50 electron masses

#### 5. **analysis/h_xe_variables.json**
- Created JSON template with all code variables and their H+Xe values

---

## Physical Parameters Used

### Toennies et al. (1979) — Designation 2 (Table V)

| Parameter | Value | SI Units | Atomic Units |
|-----------|-------|----------|--------------|
| Well Depth (ε) | 7.08 meV | – | 0.000260 Ha |
| Equilibrium Distance (r_m) | 3.82 Å | – | 7.219 a₀ |
| Switching Radius (r_s) | 4.74 Å | – | 8.957 a₀ |
| C₆ | 25.93×10³ | meV·Å⁶ | 42.30 Ha·a₀⁶ |
| C₈ | 157.73×10³ | meV·Å⁸ | 918.37 Ha·a₀⁸ |
| C₁₀ | 1325.32×10³ | meV·Å¹⁰ | 27554.64 Ha·a₀¹⁰ |

### Masses
| Atom | Mass (amu) | Electron Masses |
|------|-----------|-----------------|
| H | 1.008 | 1837.15 |
| Xe | 131.29 | 239332.50 |
| μ (reduced) | 1.0002 | 1823.16 |

---

## Expected Results

Based on Table I of Toennies et al. (1979), resonances in the H+Xe system should appear at:

| l | Energy (cm⁻¹) | Width (cm⁻¹) |
|---|---|---|
| 5 | 5.48 ± 0.32 | narrow |
| 6 | 14.52 ± 0.97 | narrow |
| 7 | 25.89 ± 1.61 | narrow |
| 8 | 39.84 ± 2.02 | narrow |

These resonances correspond to orbiting (rainbow) scattering features in the phase shifts and DCS.

---

## File Structure

```
CodesForColdCollisionProject/
├── vpa/
│   ├── Main.f90           ← Main VPA program (modified for H+Xe)
│   ├── Potential.f90      ← Potential function (H+Xe params)
│   ├── VPA-Lib.f90        ← VPA library routines
│   ├── RICBES.f           ← Riccati-Bessel functions
│   ├── DLSODALIB.f        ← ODE solver
│   ├── WKB.f90            ← WKB approximation
│   ├── cgcoef.f           ← Clebsch-Gordan coefficients
│   ├── legendre.f90       ← Legendre polynomials
│   ├── Makefile
│   ├── difelsccs          ← Compiled executable
│   ├── potential.txt      ← Output: potential grid
│   ├── phshift.txt        ← Output: phase shifts
│   └── angdcs.txt         ← Output: angular DCS
│
├── analysis/
│   ├── plot_h_xe_results.py    ← Plotting script
│   ├── h_xe_variables.json     ← Parameter reference
│   ├── variables_report.json   ← Full variable scan
│   ├── variables_report.csv    ← CSV format
│   └── plots/
│       ├── potential.png
│       ├── phshift.png
│       ├── dcs_vs_angle.png
│       └── dcs_heatmap.png
│
└── H-Kr/                  ← Reference H+Kr calculation outputs
    └── Calculation/
```

---

## Compilation Issues & Solutions

### Problem: "exec format error" or "ld: library not found for -lSystem"

**Cause**: gfortran linker using wrong system libraries (cross-compile flags).

**Solution**: Use explicit macOS SDK paths during compilation:
```bash
/opt/homebrew/bin/gfortran \
  -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk \
  -mmacosx-version-min=10.15 \
  -O2 -o difelsccs *.f *.f90
```

### Problem: Binary is ELF Linux format instead of Mach-O

**Cause**: Using pre-compiled Linux binary from repo.

**Solution**: Delete old binary and recompile:
```bash
rm -f difelsccs
# Then run compile command above
```

---

## Advanced Usage

### Different Energy Ranges

For higher-energy collisions (>50 cm⁻¹):
```bash
./difelsccs << 'EOF'
DCS
0,15,
10.0,500.0,100
0.0,180.0,180
EOF
```

### Energy Scan Only (Phase Shifts)

```bash
./difelsccs << 'EOF'
ESCAN
5
1.0,20.0
EOF
```
This calculates η_l(E) for l=5 over 1–20 cm⁻¹.

### Angular Momentum Scan

```bash
./difelsccs << 'EOF'
LSCAN
5.0
0,15
EOF
```
This calculates η_l for E=5.0 cm⁻¹ over l=0–15.

---

## References

1. Toennies, J. P., Welz, W., & Wolf, G. (1979). *J. Chem. Phys.* **71**, 614.
2. Koike, F. (2024). Hydrogen collision scattering with rare gases. Internal report.
3. Palov, A. P., & Balint-Kurti, G. G. (2020). *Comput. Phys. Commun.*, Variable Phase Approach package.

---

## Contact & Support

For questions or modifications, refer to:
- Original VPA documentation in `vpa/orig.pkgs/`
- ReadMe files in `vpa/` and repository root
- Fortran manual reference: `vpa/gfortran.manual.txt`

---

**Generated**: 2026-06-17  
**System**: H + Xe (Toennies hybrid potential)  
**Author**: Modified from original H+Kr codebase by F. Koike
