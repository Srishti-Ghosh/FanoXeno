#!/usr/bin/env python3
import subprocess
from pathlib import Path

# Paths
ANALYSIS_DIR = Path(__file__).parent
PLOTS_DIR = ANALYSIS_DIR / "plots" / "gnuplot_generated"
VPA_DIR = ANALYSIS_DIR.parent / "vpa"

POTENTIAL_FILE = VPA_DIR / "potential.txt"
WPARAM_FILE = VPA_DIR / "wparameter.txt"

def plot_effective_potentials():
    """Graph 1: Effective Potential V(r) + Centrifugal Term for l=0 to 8"""
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Reduced mass of H+Xe in electron mass units
    mu = 1823.14 
    conv = 27211.386 # Hartree to meV conversion
    
    # Build plot commands dynamically for l=0 through 8
    lines = []
    for l in range(9):
        # Formula: ( V(r) + l(l+1)/(2*mu*R^2) ) * 27211.386
        equation = f"($2 + {l}*({l}+1)/(2.0*{mu}*$1**2)) * {conv}"
        lines.append(f"'{POTENTIAL_FILE.resolve().as_posix()}' using 1:({equation}) w lines lw 2 title 'l={l}'")
    
    plot_cmd = ", \\\n         ".join(lines)
    
    script = f"""
    set terminal pdfcairo size 8,6
    set output '{PLOTS_DIR.resolve()}/koike_graph_1_effective_potentials.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Effective Interaction Potentials (l=0 to 8)"
    set xlabel "Internuclear distance R (a₀)" font "Times New Roman, 18"
    set ylabel "Effective Potential Energy (meV)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    
    # Zoom in exactly on the potential well region
    set xrange [4:20]
    set yrange [-200:20]
    
    # Draw a line at E=0 to clearly show where the well disappears
    set xzeroaxis
    
    plot {plot_cmd}
    """
    
    plt_file = PLOTS_DIR / "eff_pot.plt"
    with open(plt_file, 'w') as f: f.write(script)
    subprocess.run(['gnuplot', str(plt_file)])
    print(f"✓ Generated Graph 1: {PLOTS_DIR}/koike_graph_1_effective_potentials.pdf")


def plot_3d_w_trajectory():
    """Graph 3: 3D Angular Evolution of the Complex w(theta) Parameter with Drop Lines"""
    
    script = f"""
    set terminal pdfcairo size 8,8
    set output '{PLOTS_DIR.resolve()}/koike_graph_3_w_trajectory_3d.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: 3D Angular Evolution of w(θ) (l=7)"
    
    set xlabel "Re[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set ylabel "Im[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set zlabel "Scattering Angle θ (°)" font "Times New Roman, 16" offset 3,0,0
    
    set grid
    set border linewidth 1.5
    
    # Make the Z-axis (Angle) start exactly at the floor
    set xyplane at 0
    set zrange [0:180]
    set ztics 30
    
    # Adjust viewing angle to match the Koike 3D perspective
    set view 60, 55, 1, 1.2
    
    # High contrast color map
    set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
    
    # LAYER 1: The X-Y Plane Trace (Shadow) - Forcing Z to 0
    # LAYER 2: The Drop Lines (Impulses)
    # LAYER 3: The Thick 3D Trajectory Tube
    
    splot '{WPARAM_FILE.resolve().as_posix()}' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \\
          '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3:3 w impulses palette lw 0.5 notitle, \\
          '{WPARAM_FILE.resolve().as_posix()}' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    """
    
    plt_file = PLOTS_DIR / "w_3d.plt"
    with open(plt_file, 'w') as f: f.write(script)
    subprocess.run(['gnuplot', str(plt_file)])
    print(f"✓ Generated UPGRADED Graph 3: {PLOTS_DIR}/koike_graph_3_w_trajectory_3d.pdf")

if __name__ == "__main__":
    if not POTENTIAL_FILE.exists() or not WPARAM_FILE.exists():
        print("✗ Error: Missing Fortran data files. Run difelsccs first.")
    else:
        plot_effective_potentials()
        plot_3d_w_trajectory()