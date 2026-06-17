set term qt size 640,720
set terminal pdfcairo size 5,6
set output 'angledcs.pdf'
#set size ratio 2.00
set tics font "Times New Roman, 20"
set label font "Times New Roman, 24"
set xlabel font "Times New Roman, 24"
set ylabel font "Times New Roman, 24"
set zlabel font "Times New Roman, 24"
#set key at 5.9,2.0
set key font "Times New Roman,24"
set title font "Times New Roman, 24, bold"
#set title "Modified angle differential crosssections"
set xlabel "E" offset 0,-0.5,0
set ylabel "q" font "Symbol,24" offset 0.5,-0.5,0
set zlabel "polar-angle DCS" offset 15,17,0 
set border linewidth 1.8
set xrange [2:6] 
set yrange [0:180]
set zrange [0:5.2]
set xyplane at 0
set hidden3d
set xtics out 2,1,6 offset -0.4,-0.4,0
set ytics out 0,30,180 offset 1.8,-0.2,0
set ztics in 0,1,5
set view 63,38,1,1.25
splot 'angdcs.2pi.txt' w l palette lw 1.2 notitle, 'asydcs.2pi.txt' w l palette lw 1.2 notitle
#splot 'asydcs.raised.txt' w l lw 1.2 notitle, 'angdcs.txt' w l lw 1.2 notitle
pause -100
