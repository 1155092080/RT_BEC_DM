set term postscript eps enhanced color "Times-Roman" 30
set output "fig.eps"
set notitle
set tics out
unset key
set pm3d map
set cbtics .06 offset -1.2,0
set xtics 3 offset 0,.75
set ytics 3 offset 1,0
set mxtics 3
set mytics 3
set xrange[-6:6]
set yrange[-6:6]
set palette model RGB
set palette defined (0 0 0 0.6274, 0.25 0 1 1, 0.5 0 1 0, 0.75 1 1 0, 1 1 0 0)
set size square
set xlabel "{/Times-Italic x} ({/Symbol m}m)" offset 0,1.2
set ylabel "{{/Times-Italic y} ({/Symbol m}m)}" offset 1.2,0
set colorbox
splot 'real3d-den-2d_xy0.txt' us ($2):($1):($3)
