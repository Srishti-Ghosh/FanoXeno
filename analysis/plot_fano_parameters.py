#!/usr/bin/env python3
import subprocess
from pathlib import Path

PLOTS_DIR = Path(__file__).parent / "plots" / "gnuplot_generated"
WPARAM_FILE = Path(__file__).parent.parent / "vpa" / "wparameter.txt"
QPARAM_FILE = Path(__file__).parent.parent / "vpa" / "qparameter.txt"
ASYDCS_FILE = Path(__file__).parent.parent / "vpa" / "asydcs.txt"

def create_fano_plots():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Complex Plane w(theta) Trajectory Plot
    w_script = f"""
    set terminal pdfcairo size 6,6
    set output '{PLOTS_DIR.resolve()}/w_complex_plane.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Complex w(θ) Parameter Trajectory (l=7)"
    set xlabel "Re[w(θ)]" font "Times New Roman, 18"
    set ylabel "Im[w(θ)]" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set size square
    # Plotting Column 1 (Real) vs Column 2 (Imaginary), color-mapped to Column 3 (Angle)
    plot '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3 w lines palette lw 3 title "Angle θ Evolution"
    """
    
    w_plt = PLOTS_DIR / "w_plot.plt"
    with open(w_plt, 'w') as f: f.write(w_script)
    subprocess.run(['gnuplot', str(w_plt)])
    print(f"✓ Generated Complex w(theta) plot: {PLOTS_DIR}/w_complex_plane.pdf")

    # 2. Fano q(theta) Profile Plot
    q_script = f"""
    set terminal pdfcairo size 8,5
    set output '{PLOTS_DIR.resolve()}/q_vs_theta.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Angle-Dependent Fano q-parameter (l=7)"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "q(θ)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    plot '{QPARAM_FILE.resolve().as_posix()}' using 1:2 w lines linecolor rgb "dark-red" lw 2.5 notitle
    """
    
    q_plt = PLOTS_DIR / "q_plot.plt"
    with open(q_plt, 'w') as f: f.write(q_script)
    subprocess.run(['gnuplot', str(q_plt)])
    print(f"✓ Generated Fano q(theta) plot: {PLOTS_DIR}/q_vs_theta.pdf")

def fix_and_plot_energy_profiles():
    """Parses the data in Python first to avoid Gnuplot's NaN line-break bug."""
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    
    angles_to_extract = [30.0, 60.0, 120.0]
    extracted_data = {ang: [] for ang in angles_to_extract}

    print("📊 Filtering energy profiles from asydcs.txt...")
    # Parse file and group by our target angles
    with open(ASYDCS_FILE, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split()
                if len(parts) >= 3:
                    energy, angle, dcs = float(parts[0]), float(parts[1]), float(parts[2])
                    
                    # Float tolerance check for the target angles
                    for target in angles_to_extract:
                        if abs(angle - target) < 0.5:
                            extracted_data[target].append(f"{energy} {dcs}")

    # Write clean temporal files for Gnuplot
    plot_commands = []
    for ang in angles_to_extract:
        temp_file = PLOTS_DIR / f"slice_{int(ang)}deg.txt"
        with open(temp_file, 'w') as f:
            f.write("\n".join(extracted_data[ang]))
        plot_commands.append(f"'{temp_file.resolve().as_posix()}' using 1:2 w lines lw 2.5 title 'θ = {ang}°'")

    script = f"""
    set terminal pdfcairo size 8,6
    set output '{PLOTS_DIR.resolve()}/fano_energy_profiles.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Beutler-Fano Resonance Profiles (l=7, Er=25.89 cm⁻¹)"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    
    # Zoom tightly around the resonance
    set xrange [24.5:27.5] 
    
    plot {", ".join(plot_commands)}
    """
    
    profile_plt = PLOTS_DIR / "energy_profiles.plt"
    with open(profile_plt, 'w') as f: f.write(script)
    subprocess.run(['gnuplot', str(profile_plt)])
    print(f"✓ Generated FIXED Fano Energy Profiles: {PLOTS_DIR}/fano_energy_profiles.pdf")

def create_magnitude_and_background_plots():
    """Plots |w(theta)| and z(theta)"""
    
    # 1. Magnitude |w(theta)| Plot
    # wparameter.txt Columns: 1=Real, 2=Imag, 3=Angle
    # Magnitude = sqrt(Real^2 + Imag^2)
    w_mag_script = f"""
    set terminal pdfcairo size 8,5
    set output '{PLOTS_DIR.resolve()}/w_magnitude.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Magnitude of w-parameter |w(θ)| (l=7)"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "|w(θ)|" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    plot '{WPARAM_FILE.resolve().as_posix()}' using 3:(sqrt($1**2 + $2**2)) w lines linecolor rgb "navy" lw 2.5 notitle
    """
    w_mag_plt = PLOTS_DIR / "w_magnitude.plt"
    with open(w_mag_plt, 'w') as f: f.write(w_mag_script)
    subprocess.run(['gnuplot', str(w_mag_plt)])
    print(f"✓ Generated |w(theta)| plot: {PLOTS_DIR}/w_magnitude.pdf")

    # 2. Background Factor z(theta) Plot
    # qparameter.txt Columns: 1=Angle, 2=q(theta), 3=z(theta)
    # Log scale is usually best for cross-sections/backgrounds
    z_bg_script = f"""
    set terminal pdfcairo size 8,5
    set output '{PLOTS_DIR.resolve()}/z_background.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Background Factor z(θ) (l=7)"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "z(θ) (a.u.)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    set logscale y 
    plot '{QPARAM_FILE.resolve().as_posix()}' using 1:3 w lines linecolor rgb "forest-green" lw 2.5 notitle
    """
    z_bg_plt = PLOTS_DIR / "z_background.plt"
    with open(z_bg_plt, 'w') as f: f.write(z_bg_script)
    subprocess.run(['gnuplot', str(z_bg_plt)])
    print(f"✓ Generated z(theta) background plot: {PLOTS_DIR}/z_background.pdf")

if __name__ == "__main__":
    if not WPARAM_FILE.exists() or not QPARAM_FILE.exists():
        print("✗ Error: Run Fortran executable with 'WPRM' control word first.")
    else:
        create_fano_plots()
        fix_and_plot_energy_profiles()
        create_magnitude_and_background_plots()