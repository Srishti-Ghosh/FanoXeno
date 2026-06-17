#!/usr/bin/env python3
"""
Master Plot Generation Script for H+Xe Scattering
Generates Matplotlib and Gnuplot figures and routes them to a 'final' directory.
"""

import subprocess
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# DIRECTORY SETUP
# =============================================================================
ANALYSIS_DIR = Path(__file__).parent
VPA_DIR = ANALYSIS_DIR.parent / "vpa"
FINAL_PLOTS_DIR = ANALYSIS_DIR / "plots" / "gnuplot_generated" / "final"

# Data Files
ANGDCS_FILE = VPA_DIR / "angdcs.txt"
POTENTIAL_FILE = VPA_DIR / "potential.txt"
PHSHIFT_FILE = VPA_DIR / "phshift.txt"
WPARAM_FILE = VPA_DIR / "wparameter.txt"
QPARAM_FILE = VPA_DIR / "qparameter.txt"
ASYDCS_FILE = VPA_DIR / "asydcs.txt"
ANGDCS_2PI_FILE = FINAL_PLOTS_DIR / "angdcs.2pi.txt"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def run_gnuplot(script_content, script_name):
    """Writes a Gnuplot script to disk and executes it."""
    script_path = FINAL_PLOTS_DIR / script_name
    with open(script_path, 'w') as f:
        f.write(script_content)
    try:
        subprocess.run(['gnuplot', str(script_path)], capture_output=True, text=True, timeout=120)
    except Exception as e:
        print(f"✗ Gnuplot execution failed for {script_name}: {e}")

def convert_angdcs_to_2pi():
    """Converts angdcs.txt to .2pi format for Gnuplot 3D plotting."""
    print("📊 Preparing 2pi data format...")
    data = np.loadtxt(ANGDCS_FILE)
    unique_energies = np.unique(data[:, 0])
    
    with open(ANGDCS_2PI_FILE, 'w') as f:
        for i, energy in enumerate(unique_energies):
            mask = data[:, 0] == energy
            energy_data = data[mask]
            for angle, value in energy_data[:, [1, 2]]:
                f.write(f"{energy:.7E}   {angle:.7E}   {value:.8g}\n")
            if i < len(unique_energies) - 1:
                f.write("\n")

def load_phshift():
    data = np.loadtxt(PHSHIFT_FILE)
    return data[:, 0], data[:, 1].astype(int), data[:, 2]

def load_dcs():
    data = np.loadtxt(ANGDCS_FILE)
    return data[:, 0], data[:, 1], data[:, 2]

# =============================================================================
# MATPLOTLIB PLOTS
# =============================================================================
def plot_energy_vs_phase_shift():
    """Generates: energy_vs_phase_shift.png (formerly phshift.png)"""
    energy, l, phshift = load_phshift()
    l_values = sorted(set(l))
    
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.tab20(np.linspace(0, 1, len(l_values)))
    
    for i, l_val in enumerate(l_values):
        mask = l == l_val
        ax.plot(energy[mask], phshift[mask], 'o-', color=colors[i], label=f'l={l_val}', markersize=3, linewidth=1.5)
    
    ax.set_xlabel('Collision Energy (cm$^{-1}$)', fontsize=12)
    ax.set_ylabel('Phase Shift (radians)', fontsize=12)
    ax.set_title('H+Xe Scattering Phase Shifts', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, ncol=2, loc='best')
    
    out_file = FINAL_PLOTS_DIR / 'energy_vs_phase_shift.png'
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f"✓ Generated: {out_file.name}")

def plot_angle_vs_dcs():
    """Generates: angle_vs_dcs.png (formerly dcs_vs_angle.png)"""
    energy, angle, dcs = load_dcs()
    unique_energies = sorted(set(energy))
    selected_energies = [unique_energies[i] for i in np.linspace(0, len(unique_energies)-1, 5, dtype=int)]
    
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_energies)))
    
    for i, E_sel in enumerate(selected_energies):
        mask = np.abs(energy - E_sel) < 0.1
        if np.any(mask):
            ax.plot(angle[mask], dcs[mask], 'o-', color=colors[i], label=f'E={E_sel:.1f} cm$^{{-1}}$', markersize=4, linewidth=1.5)
    
    ax.set_xlabel('Scattering Angle (degrees)', fontsize=12)
    ax.set_ylabel('DCS (a.u.)', fontsize=12)
    ax.set_title('H+Xe Angular Differential Cross Sections', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    ax.set_xlim([0, 180])
    
    out_file = FINAL_PLOTS_DIR / 'angle_vs_dcs.png'
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f"✓ Generated: {out_file.name}")

# =============================================================================
# GNUPLOT PLOTS
# =============================================================================
def plot_gnuplot_routines():
    # 1. 3D Angular DCS (energy_angle_vs_dcs.pdf)
    script_3d = f"""
    set terminal pdfcairo size 5,6
    set output '{FINAL_PLOTS_DIR.resolve()}/energy_angle_vs_dcs.pdf'
    set tics font "Times New Roman, 20"
    set xlabel "E (cm⁻¹)" offset 0,-1 font "Times New Roman, 24"
    set ylabel "θ (degrees)" offset 1,-1 font "Times New Roman, 24"
    set zlabel "dσ/dΩ (a.u.)" offset 3,0 font "Times New Roman, 24"
    set title "H+Xe Angular DCS" font "Times New Roman, 24, bold"
    set border linewidth 1.8
    set logscale z
    set zrange [0.1:1000]
    set xrange [1:50]
    set yrange [0:180]
    set xyplane at 0.1
    set hidden3d
    set view 50, 45, 1, 1 
    splot '{ANGDCS_2PI_FILE.resolve().as_posix()}' w l palette lw 1.2 notitle
    """
    run_gnuplot(script_3d, "3d_dcs.plt")
    print(f"✓ Generated: energy_angle_vs_dcs.pdf")

    # 2. 3D Rotating GIF (energy_angle_vs_dcs_rotating.gif)
    script_gif = f"""
    set terminal gif animate delay 7 size 800,700 optimize
    set output '{FINAL_PLOTS_DIR.resolve()}/energy_angle_vs_dcs_rotating.gif'
    set tics font "Times New Roman, 14"
    set xlabel "E (cm⁻¹)" offset 0,-1 font "Times New Roman, 18"
    set ylabel "θ (degrees)" offset 1,-1 font "Times New Roman, 18"
    set zlabel "dσ/dΩ (a.u.)" offset 3,0 font "Times New Roman, 18"
    set title "H+Xe Angular DCS (Wireframe Resonance Profile)" font "Times New Roman, 20, bold"
    set border linewidth 1.5
    set logscale z
    set logscale cb
    set zrange [0.1:400]
    set cbrange [0.1:100]
    set xrange [1:50]
    set yrange [2:180]
    set xyplane at 0.1
    unset pm3d
    set hidden3d front
    set palette defined (0 "#000040", 1 "#0000ff", 2 "#00ffff", 3 "#00ff00", 4 "#ffff00", 5 "#ff0000")
    do for [ang=0:356:4] {{
        set view 60, ang, 1, 1.2
        splot '{ANGDCS_2PI_FILE.resolve().as_posix()}' every 3:3 w l palette linewidth 1.2 notitle
    }}
    set output
    """
    run_gnuplot(script_gif, "rotating_dcs.plt")
    print(f"✓ Generated: energy_angle_vs_dcs_rotating.gif")

    # 3. Effective Potentials (distance_vs_effective_potential.pdf)
    mu = 1823.14 
    conv = 27211.386 
    lines = []
    for l in range(9):
        equation = f"($2 + {l}*({l}+1)/(2.0*{mu}*$1**2)) * {conv}"
        lines.append(f"'{POTENTIAL_FILE.resolve().as_posix()}' using 1:({equation}) w lines lw 2 title 'l={l}'")
    plot_cmd = ", \\\n         ".join(lines)
    
    script_pot = f"""
    set terminal pdfcairo size 8,6
    set output '{FINAL_PLOTS_DIR.resolve()}/distance_vs_effective_potential.pdf'
    set title "H+Xe: Effective Interaction Potentials (l=0 to 8)" font "Times New Roman, 20, bold"
    set xlabel "Internuclear distance R (a₀)" font "Times New Roman, 18"
    set ylabel "Effective Potential Energy (meV)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [4:18]
    set yrange [-200:25]
    set xzeroaxis lw 3.5 dt 2 linecolor rgb "black"
    plot {plot_cmd}
    """    
    run_gnuplot(script_pot, "eff_pot.plt")
    print(f"✓ Generated: distance_vs_effective_potential.pdf")

    # 4. Fano Energy Profiles (energy_vs_dcs_profiles.pdf)
    angles_to_extract = [30.0, 60.0, 120.0]
    extracted_data = {ang: [] for ang in angles_to_extract}
    with open(ASYDCS_FILE, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split()
                if len(parts) >= 3:
                    energy, angle, dcs = float(parts[0]), float(parts[1]), float(parts[2])
                    for target in angles_to_extract:
                        if abs(angle - target) < 0.5:
                            extracted_data[target].append(f"{energy} {dcs}")

    plot_commands = []
    for ang in angles_to_extract:
        temp_file = FINAL_PLOTS_DIR / f"slice_{int(ang)}deg.txt"
        with open(temp_file, 'w') as f:
            f.write("\n".join(extracted_data[ang]))
        plot_commands.append(f"'{temp_file.resolve().as_posix()}' using 1:2 w lines lw 2.5 title 'θ = {ang}°'")

    script_fano = f"""
    set terminal pdfcairo size 8,6
    set output '{FINAL_PLOTS_DIR.resolve()}/energy_vs_dcs_profiles.pdf'
    set title "H+Xe: Beutler-Fano Resonance Profiles (l=7, Er=25.89 cm⁻¹)" font "Times New Roman, 20, bold"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [24.5:27.5] 
    plot {", ".join(plot_commands)}
    """
    run_gnuplot(script_fano, "fano_profiles.plt")
    print(f"✓ Generated: energy_vs_dcs_profiles.pdf")

    # 5. Complex Plane (real_w_vs_imag_w.pdf)
    script_w_complex = f"""
    set terminal pdfcairo size 6,6
    set output '{FINAL_PLOTS_DIR.resolve()}/real_w_vs_imag_w.pdf'
    set title "H+Xe: Complex w(θ) Parameter Trajectory (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Re[w(θ)]" font "Times New Roman, 18"
    set ylabel "Im[w(θ)]" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set size square
    plot '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3 w lines palette lw 3 title "Angle θ Evolution"
    """
    run_gnuplot(script_w_complex, "w_complex.plt")
    print(f"✓ Generated: real_w_vs_imag_w.pdf")

    # 6. Fano q(theta) Profile (angle_vs_q_parameter.pdf)
    script_q = f"""
    set terminal pdfcairo size 8,5
    set output '{FINAL_PLOTS_DIR.resolve()}/angle_vs_q_parameter.pdf'
    set title "H+Xe: Angle-Dependent Fano q-parameter (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "q(θ)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    plot '{QPARAM_FILE.resolve().as_posix()}' using 1:2 w lines linecolor rgb "dark-red" lw 2.5 notitle
    """
    run_gnuplot(script_q, "q_theta.plt")
    print(f"✓ Generated: angle_vs_q_parameter.pdf")

    # 7. Magnitude |w(theta)| (angle_vs_w_magnitude.pdf)
    script_w_mag = f"""
    set terminal pdfcairo size 8,5
    set output '{FINAL_PLOTS_DIR.resolve()}/angle_vs_w_magnitude.pdf'
    set title "H+Xe: Magnitude of w-parameter |w(θ)| (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "|w(θ)|" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    plot '{WPARAM_FILE.resolve().as_posix()}' using 3:(sqrt($1**2 + $2**2)) w lines linecolor rgb "navy" lw 2.5 notitle
    """
    run_gnuplot(script_w_mag, "w_mag.plt")
    print(f"✓ Generated: angle_vs_w_magnitude.pdf")

    # 8. Background Factor (angle_vs_z_background.pdf)
    script_z_bg = f"""
    set terminal pdfcairo size 8,5
    set output '{FINAL_PLOTS_DIR.resolve()}/angle_vs_z_background.pdf'
    set title "H+Xe: Background Factor z(θ) (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "z(θ) (a.u.)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    set logscale y 
    plot '{QPARAM_FILE.resolve().as_posix()}' using 1:3 w lines linecolor rgb "forest-green" lw 2.5 notitle
    """
    run_gnuplot(script_z_bg, "z_bg.plt")
    print(f"✓ Generated: angle_vs_z_background.pdf")

    # 9. 3D w Trajectory Tube (real_w_imag_w_vs_angle.pdf)
    script_w_3d = f"""
    set terminal pdfcairo size 8,8
    set output '{FINAL_PLOTS_DIR.resolve()}/real_w_imag_w_vs_angle.pdf'
    set title "H+Xe: 3D Angular Evolution of w(θ) (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Re[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set ylabel "Im[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set zlabel "Scattering Angle θ (°)" font "Times New Roman, 16" offset 3,0,0
    set grid
    set border linewidth 1.5
    set xyplane at 0
    set zrange [0:180]
    set ztics 30
    set view 60, 55, 1, 1.2
    set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
    splot '{WPARAM_FILE.resolve().as_posix()}' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \\
          '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3:3 w impulses palette lw 0.5 notitle, \\
          '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    """
    run_gnuplot(script_w_3d, "w_3d.plt")
    print(f"✓ Generated: real_w_imag_w_vs_angle.pdf")

    # 10. 3D w Trajectory Rotating GIF (real_w_imag_w_vs_angle_rotating.gif)
    script_wgif = f"""
    set terminal gif animate delay 5 size 800,800 optimize
    set output '{FINAL_PLOTS_DIR.resolve()}/real_w_imag_w_vs_angle_rotating.gif'
    set title "H+Xe: 3D Angular Evolution of w(θ) (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Re[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set ylabel "Im[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set zlabel "Scattering Angle θ (°)" font "Times New Roman, 16" offset 3,0,0
    set grid
    set border linewidth 1.5
    set ticslevel 0
    set xyplane at 0.1
    set zrange [0:180]
    set ztics 30
    set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
    do for [ang=0:356:4] {{
        set view 60, ang, 1, 1.2
        splot '{WPARAM_FILE.resolve().as_posix()}' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \\
          '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3:3 w impulses palette lw 0.5 notitle, \\
          '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    }}
    set output
    """
    run_gnuplot(script_wgif, "rotating_w.plt")
    print(f"✓ Generated: real_w_imag_w_vs_angle_rotating.gif")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    print("\n" + "="*70)
    print("H+Xe COLLISION SCATTERING - MASTER PLOT GENERATOR")
    print("="*70)
    
    FINAL_PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\n📁 Target Directory created/verified: {FINAL_PLOTS_DIR}")
    
    if not ANGDCS_FILE.exists() or not WPARAM_FILE.exists():
        print("✗ Error: Missing Fortran data files. Run difelsccs first.")
        sys.exit(1)
        
    convert_angdcs_to_2pi()
    
    print("\n📈 Executing Matplotlib routines...")
    plot_energy_vs_phase_shift()
    plot_angle_vs_dcs()
    
    print("\n🎨 Executing Gnuplot routines...")
    plot_gnuplot_routines()

    print("\n" + "="*70)
    print(f"✓ Success! All specific files generated in: {FINAL_PLOTS_DIR}")
    print("="*70 + "\n")