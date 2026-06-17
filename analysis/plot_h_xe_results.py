#!/usr/bin/env python3
"""
Plot H+Xe collision results from VPA code output files.
Generates figures for potential, phase shifts, and angular differential cross sections.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# Adjust paths as needed
VPA_DIR = Path(__file__).parent.parent / 'vpa'
OUTPUT_DIR = Path(__file__).parent / 'plots'
OUTPUT_DIR.mkdir(exist_ok=True)

def load_potential(filepath):
    """Load potential energy from potential.txt"""
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]  # r (Bohr), V (Hartree)

def load_phshift(filepath):
    """Load phase shifts from phshift.txt"""
    data = np.loadtxt(filepath)
    # Format: energy (cm^-1), l, phase_shift (radians)
    energy = data[:, 0]
    l = data[:, 1].astype(int)
    phshift = data[:, 2]
    return energy, l, phshift

def load_dcs(filepath):
    """Load angular differential cross sections from angdcs.txt"""
    data = np.loadtxt(filepath)
    # Format: energy (cm^-1), angle (degrees), dcs
    energy = data[:, 0]
    angle = data[:, 1]
    dcs = data[:, 2]
    return energy, angle, dcs

def plot_potential():
    """Plot the interaction potential"""
    r, V = load_potential(VPA_DIR / 'potential.txt')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(r, V * 27211.386, 'b-', linewidth=2, label='H+Xe potential')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Internuclear Distance (Bohr)', fontsize=12)
    ax.set_ylabel('Potential Energy (meV)', fontsize=12)
    ax.set_title('H+Xe Interaction Potential (Toennies et al., 1979)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    ax.set_xlim([0, 50])
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'potential.png', dpi=150)
    print(f"✓ Saved: {OUTPUT_DIR / 'potential.png'}")
    plt.close()

def plot_phshift():
    """Plot phase shifts as a function of energy for different l values"""
    energy, l, phshift = load_phshift(VPA_DIR / 'phshift.txt')
    
    # Group by l value
    l_values = sorted(set(l))
    
    fig, ax = plt.subplots(figsize=(11, 7))
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(l_values)))
    
    for i, l_val in enumerate(l_values):
        mask = l == l_val
        E_l = energy[mask]
        ps_l = phshift[mask]
        ax.plot(E_l, ps_l, 'o-', color=colors[i], label=f'l={l_val}', markersize=3, linewidth=1.5)
    
    ax.set_xlabel('Collision Energy (cm$^{-1}$)', fontsize=12)
    ax.set_ylabel('Phase Shift (radians)', fontsize=12)
    ax.set_title('H+Xe Scattering Phase Shifts', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, ncol=2, loc='best')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'phshift.png', dpi=150)
    print(f"✓ Saved: {OUTPUT_DIR / 'phshift.png'}")
    plt.close()

def plot_dcs_vs_angle_selected_energies():
    """Plot DCS vs scattering angle for selected collision energies"""
    energy, angle, dcs = load_dcs(VPA_DIR / 'angdcs.txt')
    
    # Select a subset of energies for clarity
    unique_energies = sorted(set(energy))
    selected_energies = [unique_energies[i] for i in np.linspace(0, len(unique_energies)-1, 5, dtype=int)]
    
    fig, ax = plt.subplots(figsize=(11, 7))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_energies)))
    
    for i, E_sel in enumerate(selected_energies):
        mask = np.abs(energy - E_sel) < 0.1  # small tolerance
        angle_E = angle[mask]
        dcs_E = dcs[mask]
        if len(angle_E) > 0:
            ax.plot(angle_E, dcs_E, 'o-', color=colors[i], label=f'E={E_sel:.1f} cm$^{{-1}}$', 
                   markersize=4, linewidth=1.5)
    
    ax.set_xlabel('Scattering Angle (degrees)', fontsize=12)
    ax.set_ylabel('DCS (arbitrary units)', fontsize=12)
    ax.set_title('H+Xe Angular Differential Cross Sections', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    ax.set_xlim([0, 180])
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'dcs_vs_angle.png', dpi=150)
    print(f"✓ Saved: {OUTPUT_DIR / 'dcs_vs_angle.png'}")
    plt.close()

def plot_dcs_3d_heatmap():
    """Create a 2D heatmap of DCS(energy, angle)"""
    energy, angle, dcs = load_dcs(VPA_DIR / 'angdcs.txt')
    
    # Pivot to create a 2D grid
    unique_energies = sorted(set(energy))
    unique_angles = sorted(set(angle))
    
    dcs_grid = np.zeros((len(unique_energies), len(unique_angles)))
    
    for i, E in enumerate(unique_energies):
        for j, ang in enumerate(unique_angles):
            mask = (np.abs(energy - E) < 0.01) & (np.abs(angle - ang) < 0.5)
            if np.any(mask):
                dcs_grid[i, j] = dcs[mask][0]
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    im = ax.contourf(unique_angles, unique_energies, dcs_grid, levels=50, cmap='viridis')
    cbar = plt.colorbar(im, ax=ax, label='DCS (a.u.)')
    
    ax.set_xlabel('Scattering Angle (degrees)', fontsize=12)
    ax.set_ylabel('Collision Energy (cm$^{-1}$)', fontsize=12)
    ax.set_title('H+Xe DCS: Energy vs Angle Heatmap', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'dcs_heatmap.png', dpi=150)
    print(f"✓ Saved: {OUTPUT_DIR / 'dcs_heatmap.png'}")
    plt.close()

def generate_summary_table():
    """Generate a summary table of key parameters used"""
    summary = """
H+Xe SCATTERING CALCULATION SUMMARY
====================================

Potential Model: Toennies et al. (1979) - Hybrid (L-J + 3-term dispersion)

POTENTIAL PARAMETERS (H+Xe):
  Well Depth (ε):         7.08 meV = 0.000260 Hartree
  Equilibrium Distance:   7.219 a₀ = 3.82 Å
  Switching Radius (rs):  8.957 a₀ = 4.74 Å
  
DISPERSION COEFFICIENTS:
  C₆:  1.151 × 10⁶ meV·a₀⁶  =  42.30 Hartree·a₀⁶
  C₈:  2.499 × 10⁷ meV·a₀⁸  = 918.37 Hartree·a₀⁸
  C₁₀: 7.498 × 10⁸ meV·a₀¹⁰ = 27554.64 Hartree·a₀¹⁰

MASSES:
  m_H:  1.007825 amu = 1837.15 electron masses
  m_Xe: 131.293 amu = 239332.50 electron masses
  μ_red: 1.00015 amu = 1823.16 electron masses

CALCULATION PARAMETERS:
  Angular momentum range: l = 0 to 10
  Energy range: 1.0 to 50.0 cm⁻¹
  Energy points: 200
  Angular range: 0° to 180°
  Angle points: 180

EXPECTED RESONANCE POSITIONS (from Table I, Toennies et al., 1979):
  l=5:  0.68 ± 0.04 meV  ≈ 5.48 cm⁻¹
  l=6:  1.80 ± 0.12 meV  ≈ 14.52 cm⁻¹
  l=7:  3.21 ± 0.20 meV  ≈ 25.89 cm⁻¹
  l=8:  4.94 ± 0.25 meV  ≈ 39.84 cm⁻¹

OUTPUT FILES GENERATED:
  ✓ potential.txt:  Potential energy grid
  ✓ phshift.txt:    Phase shifts vs energy and l
  ✓ angdcs.txt:     Angular DCS vs energy and angle

PLOTS GENERATED:
  ✓ potential.png:   Interaction potential curve
  ✓ phshift.png:     Phase shifts for all l values
  ✓ dcs_vs_angle.png: DCS at selected energies
  ✓ dcs_heatmap.png: 2D contour plot of DCS(E, θ)
"""
    
    summary_file = OUTPUT_DIR / 'CALCULATION_SUMMARY.txt'
    with open(summary_file, 'w') as f:
        f.write(summary)
    
    print(f"✓ Saved: {summary_file}")
    print(summary)

def main():
    print("\n" + "="*70)
    print("H+Xe COLLISION SCATTERING - PLOT GENERATION")
    print("="*70 + "\n")
    
    try:
        print("Plotting potential energy...")
        plot_potential()
        
        print("Plotting phase shifts...")
        plot_phshift()
        
        print("Plotting DCS vs angle (selected energies)...")
        plot_dcs_vs_angle_selected_energies()
        
        print("Plotting DCS heatmap...")
        plot_dcs_3d_heatmap()
        
        print("Generating summary...")
        generate_summary_table()
        
        print("\n" + "="*70)
        print(f"✓ All plots saved to: {OUTPUT_DIR}")
        print("="*70 + "\n")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print(f"Make sure you have run the VPA code to generate output files in {VPA_DIR}")
        sys.exit(1)

if __name__ == '__main__':
    main()
