#!/usr/bin/env python3
import subprocess
from pathlib import Path

PLOTS_DIR = Path(__file__).parent / "plots" / "gnuplot_generated"
WPARAM_FILE = Path(__file__).parent.parent / "vpa" / "wparameter.txt"
QPARAM_FILE = Path(__file__).parent.parent / "vpa" / "qparameter.txt"

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

if __name__ == "__main__":
    if not WPARAM_FILE.exists() or not QPARAM_FILE.exists():
        print("✗ Error: Run Fortran executable with 'WPRM' control word first.")
    else:
        create_fano_plots()