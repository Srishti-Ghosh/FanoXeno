#!/usr/bin/env python3
"""
Generate all H+Xe scattering plots (Gnuplot-based) to match H+Kr results.

This script:
1. Converts raw VPA output (angdcs.txt, potential.txt, phshift.txt) to Gnuplot-compatible formats
2. Creates Gnuplot scripts for 3D polar plots (angular DCS)
3. Runs Gnuplot to generate PDF outputs
4. Replicates all plots from the H+Kr folder

Author: Conversion script for H+Xe system
Date: 2026-06-17
"""

import subprocess
import os
import sys
from pathlib import Path
import numpy as np

# Paths
VPA_DIR = Path(__file__).parent.parent / "vpa"
ANALYSIS_DIR = Path(__file__).parent
PLOTS_DIR = ANALYSIS_DIR / "plots" / "gnuplot_generated"

def create_directories():
    """Create output directories."""
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Created plots directory: {PLOTS_DIR}")

def convert_angdcs_to_2pi(input_file, output_file):
    """
    Convert raw angdcs.txt to angdcs.2pi.txt format for Gnuplot 3D plotting.
    
    The .2pi format separates different energy sections with blank lines
    for proper Gnuplot splot command interpretation.
    """
    print(f"\n📊 Converting {input_file.name} to .2pi format...")
    
    # Read data
    data = np.loadtxt(input_file)
    
    # Get unique energies (first column)
    unique_energies = np.unique(data[:, 0])
    print(f"   Found {len(unique_energies)} energy points")
    
    # Write .2pi format with blank lines between energy sections
    with open(output_file, 'w') as f:
        for i, energy in enumerate(unique_energies):
            # Get all rows for this energy
            mask = data[:, 0] == energy
            energy_data = data[mask]
            
            # Write this energy section
            for angle, value in energy_data[:, [1, 2]]:
                f.write(f"{energy:.7E}   {angle:.7E}   {value:.8g}\n")
            
            # Add blank line between sections (except after last one)
            if i < len(unique_energies) - 1:
                f.write("\n")
    
    print(f"✓ Saved: {output_file}")
    print(f"   Energies: {unique_energies[0]:.2f} to {unique_energies[-1]:.2f} cm⁻¹")
    print(f"   Lines: {len(data)} data points")
    return len(unique_energies)

def create_gnuplot_script_angular_dcs(script_file, data_file, pdf_file, title="H+Xe Angular DCS"):
    """Create Gnuplot script for 3D angular DCS plot using log-scale safe bounds."""
    data_path = data_file.resolve().as_posix()
    pdf_path = pdf_file.resolve().as_posix()
    
    script_content = f"""set terminal pdfcairo size 5,6
set output '{pdf_path}'
set tics font "Times New Roman, 20"
set xlabel font "Times New Roman, 24"
set ylabel font "Times New Roman, 24"
set zlabel font "Times New Roman, 24"
set title font "Times New Roman, 24, bold"
set title "{title}"

# Use purely 2D character offsets to avoid log(0) coordinate crashes
set xlabel "E (cm⁻¹)" offset 0,-1
set ylabel "θ (degrees)" offset 1,-1
set zlabel "dσ/dΩ (a.u.)" offset 3,0

set border linewidth 1.8

# Apply the logarithmic scale and the zoomed bounds
set logscale z
set zrange [0.1:1000]
set xrange [1:50]
set yrange [0:180]

# Set the floor of the 3D box to match the bottom of the Z-range (NOT 0)
set xyplane at 0.1

set hidden3d
set view 50, 45, 1, 1 
splot '{data_path}' w l palette lw 1.2 notitle
"""
    with open(script_file, 'w') as f:
        f.write(script_content)
    print(f"✓ Created Gnuplot script: {script_file}")

def get_energy_range(data_file):
    """Get min and max energy from data file."""
    data = np.loadtxt(data_file)
    energies = np.unique(data[:, 0])
    return energies[0], energies[-1]

def get_dcs_max(data_file):
    """Get max DCS value for z-range."""
    data = np.loadtxt(data_file)
    return np.max(data[:, 2])

def create_gnuplot_script_potential(script_file, data_file, pdf_file, title="H+Xe Potential"):
    """Create Gnuplot script for potential energy plot using absolute paths."""
    data_path = data_file.resolve().as_posix()
    pdf_path = pdf_file.resolve().as_posix()
    
    script_content = f"""set terminal pdfcairo size 8,6
set output '{pdf_path}'
set size ratio 0.80
set tics font "Times New Roman, 18"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key font "Times New Roman,16"
set title font "Times New Roman, 24, bold"
set title "{title}"
set xlabel "Internuclear distance (a₀)"
set ylabel "Potential energy (meV)"
set border linewidth 2
set xtics 2
set ytics 2
set xrange [4:15]      # Focus on the well and immediate long-range tail
set yrange [-10:5]     # Clip the repulsive wall to see the -7.08 meV well
plot '{data_path}' using 1:(($2)*27211.386) w l lt 8 lw 3 title "V(r) H+Xe"
"""
    with open(script_file, 'w') as f:
        f.write(script_content)
    print(f"✓ Created Gnuplot script: {script_file}")

def create_gnuplot_script_phshift(script_file, data_file, pdf_file, title="H+Xe Phase Shifts"):
    """Create Gnuplot script for phase shift plot dynamically handling all partial waves."""
    data_path = data_file.resolve().as_posix()
    pdf_path = pdf_file.resolve().as_posix()
    
    data = np.loadtxt(data_file)
    num_cols = data.shape[1]
    
    plot_lines = []
    for col in range(2, num_cols + 1):
        l_val = col - 2
        plot_lines.append(f"'{data_path}' using 1:{col} w l lw 2.2 title 'η_{l_val}(E)'")
    plot_cmd = "plot " + ", \\\n     ".join(plot_lines)

    script_content = f"""set terminal pdfcairo size 8,6
set output '{pdf_path}'
set tics font "Times New Roman, 16"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key outside right top font "Times New Roman,14"
set title font "Times New Roman, 24, bold"
set title "{title}"
set xlabel "Energy (cm⁻¹)"
set ylabel "Phase Shift η_l (radians)"
set border linewidth 2
set xtics 5
set xrange [0:50]
set yrange [-5:10]      # Adjust based on your data's specific phase range
set ytics 3.14159       # Set tick marks to multiples of Pi (π)
set format y "%.1f"
{plot_cmd}
"""
    with open(script_file, 'w') as f:
        f.write(script_content)
    print(f"✓ Created Gnuplot script: {script_file}")

def run_gnuplot(script_file):
    """Run Gnuplot script to generate PDF."""
    try:
        result = subprocess.run(['gnuplot', str(script_file)], 
                              capture_output=True, 
                              text=True,
                              timeout=60) # Increased timeout for font caching safety
        if result.returncode == 0:
            pdf_name = script_file.name.replace('.plt', '.pdf')
            print(f"   ✓ Generated: {pdf_name}")
            return True
        else:
            print(f"   ✗ Gnuplot error:\n{result.stderr}")
            return False
    except FileNotFoundError:
        print("   ✗ Gnuplot binary execution failed.")
        return False
    except subprocess.TimeoutExpired:
        print("   ✗ Gnuplot processing timed out")
        return False

def check_gnuplot_installed():
    """Check if Gnuplot is installed with an extended timeout safety margin."""
    try:
        subprocess.run(['gnuplot', '--version'], 
                      capture_output=True, 
                      timeout=30)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False

def main():
    print("=" * 70)
    print("H+XE GNUPLOT-BASED PLOT GENERATION")
    print("=" * 70)
    
    create_directories()
    
    if not check_gnuplot_installed():
        print("\n⚠️  WARNING: Gnuplot execution timed out or not found!")
        print("Please check your local installation or system environment path variables.")
        gnuplot_available = False
    else:
        print("✓ Gnuplot is installed and responding")
        gnuplot_available = True
    
    angdcs_file = VPA_DIR / "angdcs.txt"
    potential_file = VPA_DIR / "potential.txt"
    phshift_file = VPA_DIR / "phshift.txt"
    
    if not angdcs_file.exists():
        print(f"✗ Error: {angdcs_file} not found. Run your VPA calculation engine first.")
        return 1
    
    angdcs_2pi = PLOTS_DIR / "angdcs.2pi.txt"
    convert_angdcs_to_2pi(angdcs_file, angdcs_2pi)
    
    if gnuplot_available:
        print("\n📈 Creating Gnuplot scripts...")
        
        # 1. Angular DCS
        gp_angdcs = PLOTS_DIR / "angledcs.plt"
        pdf_angdcs = PLOTS_DIR / "angledcs.pdf"
        create_gnuplot_script_angular_dcs(gp_angdcs, angdcs_2pi, pdf_angdcs)
        
        print(f"\n🎨 Generating PDFs with Gnuplot...")
        run_gnuplot(gp_angdcs)
        
        # 2. Potential Plot
        if potential_file.exists():
            gp_pot = PLOTS_DIR / "pot.plt"
            pdf_pot = PLOTS_DIR / "pot.pdf"
            create_gnuplot_script_potential(gp_pot, potential_file, pdf_pot)
            run_gnuplot(gp_pot)
        else:
            print(f"⚠️  Skipping potential graph; {potential_file.name} is missing.")
        
        # 3. Phase Shift Plot
        if phshift_file.exists():
            gp_ps = PLOTS_DIR / "psplot.plt"
            pdf_ps = PLOTS_DIR / "psplot.pdf"
            create_gnuplot_script_phshift(gp_ps, phshift_file, pdf_ps)
            run_gnuplot(gp_ps)
        else:
            print(f"⚠️  Skipping phase shift graph; {phshift_file.name} is missing.")
    
    print("\n" + "=" * 70)
    print("✓ Plot generation complete!")
    print(f"   Data and script formats saved to: {PLOTS_DIR}")
    print("=" * 70)
    return 0

if __name__ == "__main__":
    sys.exit(main())