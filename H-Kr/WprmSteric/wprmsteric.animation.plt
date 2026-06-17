set terminal gif animate delay 1 size 640,720
set output "wpmanimation.gif"
#set term qt size 640,720
#set terminal pdfcairo
#set output 'angledcs.pdf'
#set size ratio 2.00
set tics font "Times New Roman, 18"
set label font "Times New Roman, 24"
set xlabel font "Times New Roman, 24"
set ylabel font "Times New Roman, 24"
set zlabel font "Times New Roman, 24"
#set key at 5.9,2.0
set key font "Times New Roman,24"
set title font "Times New Roman, 24, bold"
set xlabel "x" font "Times New Roman, 32" offset 1.5,-0.5,0 
set ylabel "y" font "Times New Roman, 32" offset -1.5,-0.5,0
#set zlabel "q" font "Symbol, 36" offset 12,13,0
set border linewidth 1.8
set xrange [-20:40] 
set yrange [-40:20]
set zrange [0:180]
set xyplane at 0
set hidden3d
set xtics in -10,10,40 offset  1.4,-0.4,0
set ytics in -40,10,10 offset -1.8,-0.2,0
set ztics in 0,30,180
do for [i=0:360] {
set view 66, i, 1, 1.25
splot 'wparameter.txt' w l palette lw 6 notitle, 'wprmaux.txt' w l lw 4 lc 2 notitle, \
      'wparameter.txt' w i palette lw 0.5 notitle, \
      'xaxis.txt' w l lw 1.8 lc 8 notitle, 'yaxis.txt' w l lw 1.8 lc 8 notitle
#pause 0.05
}
#pause -100
set output

