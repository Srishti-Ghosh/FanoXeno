set terminal pdfcairo size 5,6
set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/angledcs.pdf'
set tics font "Times New Roman, 20"
set xlabel font "Times New Roman, 24"
set ylabel font "Times New Roman, 24"
set zlabel font "Times New Roman, 24"
set title font "Times New Roman, 24, bold"
set title "H+Xe Angular DCS"

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
splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/angdcs.2pi.txt' w l palette lw 1.2 notitle
