#!/usr/bin/env python3
"""
Master Automation & Plot Generation Script for H+Xe Scattering
Separates Global Macroscopic Physics from Local Fano Resonance Physics.
Includes dynamic auto-scaling and l-dependent angle slicing.
"""

import subprocess
import sys
import os
import shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# DIRECTORY SETUP
# =============================================================================
ANALYSIS_DIR = Path(__file__).parent
VPA_DIR = ANALYSIS_DIR.parent / "vpa"
GLOBAL_DIR = ANALYSIS_DIR / "Global_System_Plots"

# =============================================================================
# 1. FORTRAN AUTOMATION & PARSING
# =============================================================================
def parse_input_and_run(input_filepath):
    with open(input_filepath, 'r') as f:
        raw_lines = f.readlines()
        
    lines = [line.split('!')[0].strip() for line in raw_lines if line.split('!')[0].strip()]
    
    l_params = lines[1].split(',')
    lmin, lmax, lr = int(l_params[0]), int(l_params[1]), int(l_params[2])
    
    e_params = lines[2].split(',')
    emin_glob = float(e_params[0].lower().replace('d', 'e'))
    emax_glob = float(e_params[1].lower().replace('d', 'e'))
    nestep_glob = int(e_params[2])
    
    er_params = lines[3].split(',')
    er = float(er_params[0].lower().replace('d', 'e'))
    width = float(er_params[1].lower().replace('d', 'e'))
    
    ang_params = lines[4].split(',')
    angmin, angmax, nangstep = ang_params[0].strip(), ang_params[1].strip(), ang_params[2].strip()
    
    ps_block = "\n".join(lines[5:])
    
    LOCAL_DIR = ANALYSIS_DIR / f"l_{lr}_plots"
    GLOBAL_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n⚙️ Preparing Fortran engine...")
    mac_fix = 'export LIBRARY_PATH="$LIBRARY_PATH:$(xcrun --show-sdk-path)/usr/lib" && '
    subprocess.run(mac_fix + 'make', shell=True, cwd=VPA_DIR, capture_output=True)

    # --- RUN 1: GLOBAL MACROSCOPIC PHYSICS ---
    print("   🌐 Running Global Physics (0 to 50 cm⁻¹)...")
    temp_global = VPA_DIR / "temp_global.txt"
    with open(temp_global, 'w') as f:
        f.write(f"DCS\n{lmin}, {lmax}\n{emin_glob}d0, {emax_glob}d0, {nestep_glob}\n{angmin}, {angmax}, {nangstep}\n")
    
    with open(temp_global, 'r') as stdin_file:
        subprocess.run(['./difelsccs'], stdin=stdin_file, cwd=VPA_DIR)
        
    shutil.move(VPA_DIR / "potential.txt", GLOBAL_DIR / "potential.txt")
    shutil.move(VPA_DIR / "phshift.txt", GLOBAL_DIR / "phshift.txt")
    shutil.move(VPA_DIR / "angdcs.txt", GLOBAL_DIR / "angdcs.txt")

    # --- RUN 2: LOCAL FANO PHYSICS ---
    print(f"   🔬 Running Fano Resonance Physics (Zoomed at E={er} cm⁻¹)...")
    temp_wprm = VPA_DIR / "temp_wprm.txt"
    zoom_emin, zoom_emax = er - 2.0, er + 2.0
    with open(temp_wprm, 'w') as f:
        f.write(f"WPRM\n{lmin}, {lmax}, {lr}\n{zoom_emin}d0, {zoom_emax}d0, 200\n")
        f.write(f"{er}d0, {width}d0\n{angmin}, {angmax}, {nangstep}\n{ps_block}\n")

    with open(temp_wprm, 'r') as stdin_file:
        subprocess.run(['./difelsccs'], stdin=stdin_file, cwd=VPA_DIR)

    shutil.move(VPA_DIR / "wparameter.txt", LOCAL_DIR / "wparameter.txt")
    shutil.move(VPA_DIR / "qparameter.txt", LOCAL_DIR / "qparameter.txt")
    shutil.move(VPA_DIR / "asydcs.txt", LOCAL_DIR / "asydcs.txt")

    if temp_global.exists(): os.remove(temp_global)
    if temp_wprm.exists(): os.remove(temp_wprm)
        
    print("✓ Fortran calculations complete.")
    return lr, er, emin_glob, emax_glob, zoom_emin, zoom_emax, LOCAL_DIR

def convert_to_2pi(source_file, target_file):
    data = np.loadtxt(source_file)
    unique_energies = np.unique(data[:, 0])
    with open(target_file, 'w') as f:
        for i, energy in enumerate(unique_energies):
            mask = data[:, 0] == energy
            for angle, value in data[mask][:, [1, 2]]:
                f.write(f"{energy:.7E}   {angle:.7E}   {value:.8g}\n")
            if i < len(unique_energies) - 1:
                f.write("\n")

def run_gnuplot(script_content, script_name, target_dir):
    script_path = target_dir / script_name
    with open(script_path, 'w') as f: f.write(script_content)
    subprocess.run(['gnuplot', str(script_path)], capture_output=True)

# =============================================================================
# 2. PLOT GLOBAL MACROSCOPIC DATA
# =============================================================================
def plot_global_data(emin, emax):
    print("\n📈 Generating Global System Plots...")
    
    # 1. Phase Shifts (Matplotlib)
    data = np.loadtxt(GLOBAL_DIR / "phshift.txt")
    energy, l, phshift = data[:, 0], data[:, 1].astype(int), data[:, 2]
    l_values = sorted(set(l))
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.tab20(np.linspace(0, 1, len(l_values)))
    for i, l_val in enumerate(l_values):
        mask = l == l_val
        ax.plot(energy[mask], phshift[mask], 'o-', color=colors[i], label=f'l={l_val}', markersize=3, linewidth=1.5)
    ax.set_xlabel('Collision Energy (cm$^{-1}$)', fontsize=14)
    ax.set_ylabel('Phase Shift (radians)', fontsize=14)
    ax.set_title('H+Xe Scattering Phase Shifts', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10, ncol=2)
    plt.tight_layout(); plt.savefig(GLOBAL_DIR / 'energy_vs_phase_shift.png', dpi=150); plt.close()

    # 2. Angle vs Standard DCS (Matplotlib)
    data = np.loadtxt(GLOBAL_DIR / "angdcs.txt")
    energy, angle, dcs = data[:, 0], data[:, 1], data[:, 2]
    unique_energies = sorted(set(energy))
    selected_energies = [unique_energies[i] for i in np.linspace(0, len(unique_energies)-1, 5, dtype=int)]
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_energies)))
    for i, E_sel in enumerate(selected_energies):
        mask = np.abs(energy - E_sel) < 0.1
        if np.any(mask):
            ax.plot(angle[mask], dcs[mask], 'o-', color=colors[i], label=f'E={E_sel:.1f} cm$^{{-1}}$', markersize=4)
    ax.set_xlabel('Scattering Angle (degrees)', fontsize=14)
    ax.set_ylabel('DCS (a.u.)', fontsize=14)
    ax.set_title('H+Xe Standard Angular Differential Cross Sections', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3); ax.legend(fontsize=10); ax.set_xlim([0, 180])
    plt.tight_layout(); plt.savefig(GLOBAL_DIR / 'angle_vs_dcs.png', dpi=150); plt.close()

    # 3 & 4. 3D Standard DCS and GIF (Gnuplot)
    convert_to_2pi(GLOBAL_DIR / "angdcs.txt", GLOBAL_DIR / "angdcs.2pi.txt")
    common_3d = f"""
    set tics font "Times New Roman, 16"
    set xlabel "E (cm⁻¹)" offset 0,-1 font "Times New Roman, 20"
    set ylabel "θ (degrees)" offset 1,-1 font "Times New Roman, 20"
    set zlabel "dσ/dΩ (a.u.)" offset 3,0 font "Times New Roman, 20"
    set border linewidth 1.5; set logscale z; set zrange [0.1:1000]
    set xrange [{emin}:{emax}]; set yrange [0:180]; set xyplane at 0.1
    set palette defined (0 "#000040", 1 "#0000ff", 2 "#00ffff", 3 "#00ff00", 4 "#ffff00", 5 "#ff0000")
    """
    script_3d = f"set terminal pdfcairo size 6,6\nset output '{GLOBAL_DIR}/energy_angle_vs_dcs.pdf'\n" + common_3d + f"""
    set title "H+Xe Standard Angular DCS" font "Times New Roman, 20, bold"
    set hidden3d front; set view 50, 45, 1, 1 
    splot '{GLOBAL_DIR}/angdcs.2pi.txt' w l palette lw 1.2 notitle
    """
    run_gnuplot(script_3d, "gp_3d_dcs.plt", GLOBAL_DIR)

    script_gif = f"set terminal gif animate delay 7 size 800,700 optimize\nset output '{GLOBAL_DIR}/energy_angle_vs_dcs_rotating.gif'\n" + common_3d + f"""
    set title "H+Xe Standard Angular DCS" font "Times New Roman, 20, bold"
    unset pm3d; set hidden3d front; set logscale cb; set cbrange [0.1:100]
    do for [ang=0:356:4] {{ set view 60, ang, 1, 1.2; splot '{GLOBAL_DIR}/angdcs.2pi.txt' every 3:3 w l palette lw 1.2 notitle }}
    """
    run_gnuplot(script_gif, "gp_dcs_gif.plt", GLOBAL_DIR)

    # 5. Effective Potentials
    lines = [f"'{GLOBAL_DIR}/potential.txt' using 1:(($2 + {l}*({l}+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l={l}'" for l in range(11)]
    script_pot = f"""
    set terminal pdfcairo size 8,6; set output '{GLOBAL_DIR}/distance_vs_effective_potential.pdf'
    set title "H+Xe: Effective Interaction Potentials" font "Times New Roman, 20, bold"
    set xlabel "Internuclear distance R (a₀)" font "Times New Roman, 18"
    set ylabel "Effective Potential Energy (meV)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [4:18]; set yrange [-200:20]
    set xzeroaxis lw 1.5 dt 2 linecolor rgb "black"
    plot {", ".join(lines)}
    """    
    run_gnuplot(script_pot, "gp_eff_pot.plt", GLOBAL_DIR)

# =============================================================================
# 3. PLOT LOCAL FANO DATA
# =============================================================================
def plot_local_data(LOCAL_DIR, lr, er, z_emin, z_emax):
    print(f"🎨 Executing Gnuplot routines for Fano l={lr}...")
    
    asydcs_2pi = LOCAL_DIR / "asydcs.2pi.txt"
    convert_to_2pi(LOCAL_DIR / "asydcs.txt", asydcs_2pi)
    w_file = (LOCAL_DIR / "wparameter.txt").resolve().as_posix()
    q_file = (LOCAL_DIR / "qparameter.txt").resolve().as_posix()

    common_3d = f"""
    set tics font "Times New Roman, 16"
    set xlabel "E (cm⁻¹)" offset 0,-1 font "Times New Roman, 20"
    set ylabel "θ (degrees)" offset 1,-1 font "Times New Roman, 20"
    set zlabel "dσ/dΩ (a.u.)" offset 3,0 font "Times New Roman, 20"
    set border linewidth 1.5; set logscale z; set zrange [0.1:1000]
    set xrange [{z_emin}:{z_emax}]; set yrange [0:180]; set xyplane at 0.1
    set palette defined (0 "#000040", 1 "#0000ff", 2 "#00ffff", 3 "#00ff00", 4 "#ffff00", 5 "#ff0000")
    """

    # 1. & 2. 3D Fano DCS & GIF
    script_3d = f"set terminal pdfcairo size 6,6\nset output '{LOCAL_DIR}/energy_angle_vs_asydcs.pdf'\n" + common_3d + f"""
    set title "H+Xe Synthesized Fano DCS (l={lr})" font "Times New Roman, 20, bold"
    set hidden3d front; set view 50, 45, 1, 1 
    splot '{asydcs_2pi.resolve().as_posix()}' w l palette lw 1.2 notitle
    """
    run_gnuplot(script_3d, "gp_3d_asydcs.plt", LOCAL_DIR)

    script_gif = f"set terminal gif animate delay 7 size 800,700 optimize\nset output '{LOCAL_DIR}/energy_angle_vs_asydcs_rotating.gif'\n" + common_3d + f"""
    set title "H+Xe Synthesized Fano DCS (l={lr})" font "Times New Roman, 20, bold"
    unset pm3d; set hidden3d front; set logscale cb; set cbrange [0.1:100]
    do for [ang=0:356:4] {{ set view 60, ang, 1, 1.2; splot '{asydcs_2pi.resolve().as_posix()}' every 3:3 w l palette lw 1.2 notitle }}
    """
    run_gnuplot(script_gif, "gp_asydcs_gif.plt", LOCAL_DIR)

    # 3. Fano Energy Profiles (DYNAMICALLY SLICED BY FANO ARCHETYPES)
    # Load q-parameters to mathematically hunt the 4 archetypal Fano shapes
    q_data = np.loadtxt(q_file)
    angles_q, q_vals = q_data[:, 0], q_data[:, 1]
    idx_peak = np.argmax(np.abs(q_vals))    # Archetype 1: Pure Peak (q approaches infinity / maximum absolute value)
    idx_dip = np.argmin(np.abs(q_vals))    # Archetype 2: Pure Window/Dip (q approaches 0)
    idx_pos = np.argmin(np.abs(q_vals - 1.5))    # Archetype 3: Positive Asymmetry (q ≈ +1.5)
    idx_neg = np.argmin(np.abs(q_vals + 1.5))    # Archetype 4: Negative Asymmetry (q ≈ -1.5)
    target_angles = {angles_q[idx_peak], angles_q[idx_dip], angles_q[idx_pos], angles_q[idx_neg]}    # Extract these specific target angles
    extracted_data = {ang: [] for ang in target_angles}    
    with open(LOCAL_DIR / "asydcs.txt", 'r') as f:
        for line in f:
            if line.strip():
                p = line.split()
                if len(p) >= 3:
                    file_ang = float(p[1])
                    for t_ang in extracted_data.keys():
                        # Match with a tight 0.5-degree tolerance
                        if abs(file_ang - t_ang) < 0.5: 
                            extracted_data[t_ang].append(f"{p[0]} {p[2]}")
    cmds = []
    for ang in sorted(extracted_data.keys()):
        if extracted_data[ang]:  # Ensure data was actually found
            q_label = q_vals[np.argmin(np.abs(angles_q - ang))]
            tmp = LOCAL_DIR / f"slice_{int(ang)}.txt"
            with open(tmp, 'w') as f: 
                f.write("\n".join(extracted_data[ang]))
            cmds.append(f"'{tmp.resolve().as_posix()}' using 1:2 w lines lw 2.5 title 'θ = {int(ang)}° (q ≈ {q_label:.1f})'")
    script_fano = f"""
    set terminal pdfcairo size 8,6; set output '{LOCAL_DIR}/energy_vs_dcs_profiles.pdf'
    set title "H+Xe: Archetypal Fano Resonance Profiles (l={lr}, Er={er} cm⁻¹)" font "Times New Roman, 20, bold"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [{er - 1.5}:{er + 1.5}] 
    plot {", ".join(cmds)}
    """
    run_gnuplot(script_fano, "gp_fano_prof.plt", LOCAL_DIR)

    # 4. Complex Plane 2D
    script_w_comp = f"""
    set terminal pdfcairo size 6,6; set output '{LOCAL_DIR}/real_w_vs_imag_w.pdf'
    set title "H+Xe: Complex w(θ) Trajectory (l={lr})" font "Times New Roman, 20, bold"
    set xlabel "Re[w(θ)]" font "Times New Roman, 18"; set ylabel "Im[w(θ)]" font "Times New Roman, 18"
    set grid; set border linewidth 2; set size square
    plot '{w_file}' using 1:2:3 w lines palette lw 3 title "Angle θ"
    """
    run_gnuplot(script_w_comp, "gp_w_comp.plt", LOCAL_DIR)

    # 5. Fano q(theta)
    script_q = f"""
    set terminal pdfcairo size 8,5; set output '{LOCAL_DIR}/angle_vs_q_parameter.pdf'
    set title "H+Xe: Angle-Dependent Fano q-parameter (l={lr})" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"; set ylabel "q(θ)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [0:180]; set xtics 30
    plot '{q_file}' using 1:2 w lines linecolor rgb "dark-red" lw 2.5 notitle
    """
    run_gnuplot(script_q, "gp_q_theta.plt", LOCAL_DIR)

    # 6. Magnitude |w(theta)|
    script_w_mag = f"""
    set terminal pdfcairo size 8,5; set output '{LOCAL_DIR}/angle_vs_w_magnitude.pdf'
    set title "H+Xe: Magnitude |w(θ)| (l={lr})" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"; set ylabel "|w(θ)|" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [0:180]; set xtics 30
    plot '{w_file}' using 3:(sqrt($1**2 + $2**2)) w lines linecolor rgb "navy" lw 2.5 notitle
    """
    run_gnuplot(script_w_mag, "gp_w_mag.plt", LOCAL_DIR)

    # 7. Background Factor z(theta)
    script_z_bg = f"""
    set terminal pdfcairo size 8,5; set output '{LOCAL_DIR}/angle_vs_z_background.pdf'
    set title "H+Xe: Background Factor z(θ) (l={lr})" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"; set ylabel "z(θ) (a.u.)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [0:180]; set xtics 30; set logscale y 
    plot '{q_file}' using 1:3 w lines linecolor rgb "forest-green" lw 2.5 notitle
    """
    run_gnuplot(script_z_bg, "gp_z_bg.plt", LOCAL_DIR)

    # DYNAMIC BOX CALCULATION FOR 3D W-TRAJECTORY
    real_w, imag_w = [], []
    with open(w_file, 'r') as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                try:
                    real_w.append(float(p[0]))
                    imag_w.append(float(p[1]))
                except: pass
    min_x, max_x = min(real_w), max(real_w)
    min_y, max_y = min(imag_w), max(imag_w)
    
    # Calculate a perfect, tight cube base for the animation
    span = max((max_x - min_x), (max_y - min_y)) * 1.05
    cx, cy = (max_x + min_x)/2, (max_y + min_y)/2
    x_sq_min, x_sq_max = cx - span/2, cx + span/2
    y_sq_min, y_sq_max = cy - span/2, cy + span/2

    w_3d_common = f"""
    set title "H+Xe: 3D Angular Evolution of w(θ) (l={lr})" font "Times New Roman, 20, bold"    
    set xlabel "Re[w(θ)]" font "Times New Roman, 16" offset 0,-2,0
    set ylabel "Im[w(θ)]" font "Times New Roman, 16" offset 0,-2,0
    set zlabel "Scattering Angle θ (°)" font "Times New Roman, 16" offset 3,0,0
    set grid; set border linewidth 1.5; set xyplane at 0; set zrange [0:180]; set ztics 30    
    set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
    """
    
    # 8. 3D w Trajectory (Static PDF) - Autoscaled naturally!
    script_w_3d = f"set terminal pdfcairo size 8,8\nset output '{LOCAL_DIR}/real_w_imag_w_vs_angle.pdf'\n" + w_3d_common + f"""
    set view 60, 55, 1, 1.2
    splot '{w_file}' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \\
          '{w_file}' using 1:2:3:3 w impulses palette lw 0.5 notitle, \\
          '{w_file}' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    """
    run_gnuplot(script_w_3d, "gp_w_3d.plt", LOCAL_DIR)

    # 9. 3D w Trajectory (Rotating GIF) - Rigid Perfect Square Base to prevent wobbling!
    script_w_anim = f"set terminal gif animate delay 7 size 800,800 optimize\nset output '{LOCAL_DIR}/real_w_imag_w_vs_angle_rotating.gif'\n" + w_3d_common + f"""
    set xrange [{x_sq_min}:{x_sq_max}]
    set yrange [{y_sq_min}:{y_sq_max}]
    do for [ang=0:356:4] {{
        set view 60, ang, 1, 1.2
        splot '{w_file}' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \\
              '{w_file}' using 1:2:3:3 w impulses palette lw 0.5 notitle, \\
              '{w_file}' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    }}
    """
    run_gnuplot(script_w_anim, "gp_w_anim.plt", LOCAL_DIR)

# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_final_plots.py <path_to_input_file>")
        sys.exit(1)

    input_filepath = Path(sys.argv[1])
    if not input_filepath.exists():
        print(f"✗ Error: Input file '{input_filepath}' not found.")
        sys.exit(1)

    print("\n" + "="*70)
    print(f"H+Xe COLLISION SCATTERING - MASTER DUAL-RUN GENERATOR")
    print("="*70)
    
    lr, er, emin, emax, z_emin, z_emax, local_dir = parse_input_and_run(input_filepath)
    
    plot_global_data(emin, emax)
    plot_local_data(local_dir, lr, er, z_emin, z_emax)

    print("\n" + "="*70)
    print(f"✓ Success!")
    print(f"   Macroscopic System Data -> {GLOBAL_DIR.name}")
    print(f"   Fano Resonance Data     -> {local_dir.name}")
    print("="*70 + "\n")