#!/usr/bin/env python3
"""
Generate a smooth rotating GIF animation of the 3D H+Xe Angular DCS surface.
Updated with pure colored wireframe (Hidden3D) to exactly match Koike et al.
"""

import subprocess
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent
PLOTS_DIR = ANALYSIS_DIR / "plots" / "gnuplot_generated"
DATA_FILE = PLOTS_DIR / "angdcs.2pi.txt"
GIF_FILE = PLOTS_DIR / "angledcs_rotating.gif"

def create_gif_engine():
    if not DATA_FILE.exists():
        print(f"✗ Error: Complete {DATA_FILE.name} first before running animation.")
        return False

    print("🎬 Constructing Gnuplot Wireframe rotation routine...")
    data_path = DATA_FILE.resolve().as_posix()
    gif_path = GIF_FILE.resolve().as_posix()

    gnuplot_commands = f"""
    # Slowed down the animation: delay 7 instead of 5
    set terminal gif animate delay 7 size 800,700 optimize
    set output '{gif_path}'
    
    set tics font "Times New Roman, 14"
    set xlabel font "Times New Roman, 18"
    set ylabel font "Times New Roman, 18"
    set zlabel font "Times New Roman, 18"
    set title font "Times New Roman, 20, bold"
    set title "H+Xe Angular DCS (Wireframe Resonance Profile)"
    
    set xlabel "E (cm⁻¹)" offset 0,-1
    set ylabel "θ (degrees)" offset 1,-1
    set zlabel "dσ/dΩ (a.u.)" offset 3,0
    
    set border linewidth 1.5
    
    set logscale z
    set logscale cb
    set zrange [0.1:400]
    set cbrange [0.1:100]
    set xrange [1:50]
    set yrange [2:180]
    set xyplane at 0.1
    
    # THE FIX: Pure Wireframe with Hidden Lines removed
    unset pm3d
    set hidden3d front
    
    # High-contrast heat map palette
    set palette defined (0 "#000040", 1 "#0000ff", 2 "#00ffff", 3 "#00ff00", 4 "#ffff00", 5 "#ff0000")
    
    do for [ang=0:356:4] {{
        set view 60, ang, 1, 1.2
        
        # 'every 3:3' skips data points to make the grid lines clearly visible!
        # 'w l palette' colors the wireframe based on the Z-axis height
        splot '{data_path}' every 3:3 w l palette linewidth 1.2 notitle
    }}
    set output
    """
    
    plt_script = PLOTS_DIR / "animate_dcs.plt"
    with open(plt_script, 'w') as f:
        f.write(gnuplot_commands)
        
    print("🎨 Rendering Wireframe frames into GIF (this will take about 30 seconds)...")
    try:
        result = subprocess.run(['gnuplot', str(plt_script)], capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print(f"✓ Success! Animated file generated at:\n   {GIF_FILE.resolve().as_posix()}")
            return True
        else:
            print(f"✗ Gnuplot Animation Error:\n{result.stderr}")
            return False
    except Exception as e:
        print(f"✗ Failed to run process: {e}")
        return False

if __name__ == "__main__":
    create_gif_engine()