
    # Slowed down the animation: delay 7 instead of 5
    set terminal gif animate delay 7 size 800,700 optimize
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/angledcs_rotating.gif'
    
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
    
    do for [ang=0:356:4] {
        set view 60, ang, 1, 1.2
        
        # 'every 3:3' skips data points to make the grid lines clearly visible!
        # 'w l palette' colors the wireframe based on the Z-axis height
        splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/angdcs.2pi.txt' every 3:3 w l palette linewidth 1.2 notitle
    }
    set output
    