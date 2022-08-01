set terminal postscript colour font "Times-Roman,15"
set output "figure_8.ps"
set font "Times-Roman,30"
set multiplot
set label 'a)' at screen 0.18,0.90
set label 'b)' at screen 0.18,0.45
set label 'c)' at screen 0.54,0.45

set origin 0.17,0.0
set size 0.38,0.5  
unset xrange

set ylabel "{/Symbol-Oblique r}_2"
set xlabel "{/Times-Italic t}/fs"
set xrange [0:500]
set key bottom right
plot "fig_8_density_matrix_03.txt" u ($1*0.024188):2 with lines lw 7 lc rgb "#9e0142" title "{/Symbol-Oblique l}=0.3{/Symbol-Oblique l}_0", "fig_8_density_matrix_13.txt" u ($1*0.024188):2 with lines lw 7 lc rgb "#3288bd" title "{/Symbol-Oblique l}=1.3{/Symbol-Oblique l}_0"

set origin 0.53,-0.03
set size 0.39,0.53  

set ytics 0.15,0.1,0.65
set xtics 0,0.4,1.6
set xrange [0.0:1.6]
set ylabel "{/Symbol-Oblique r}_R({/Times-Italic t}={/Times-Italic t_f})"
set xlabel "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0"
plot "fig_8_yield.txt" u 1:(1-$6) with lines lw 4 lc "black" notitle, "fig_8_yield.txt" u 1:(1-$6) ps 1 pt 5 lc "black" notitle

set origin 0.20,0.50
set size 0.76,0.50

set ytics -1,0.5,1
set xtics -2,0.5,2
unset colorbox
# Conical intersection parameters
Q0 = (4.2538666725158691 + 4.8292193412780762)/2.0
kappa1 = -0.003858675*3.0
kappa2 = 0.0054756457*2.4
omegasc = 0.0043364174
omegast = 0.0027194482
e1 = 0.144792242
e2 = 0.177866612*0.87
lambda = 0.0096283166*1.3
a = 0.123646
# Lower and upper adiabatic surfaces of the conical intersection
f(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 - sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))

set xlabel "{/Times-Italic Q}_t/{/Times-Italic Q}_0"
set xtics
set ytics
#set palette defined (10 '#ffffff', 20 '#d1e5f0', 30 '#92c5de', 40 '#4393c3', 50 '#2166ac', 60 '#053061') # for nonlog
set palette defined (20 '#053061', 30 '#2166ac', 40 '#4393c3', 50 '#92c5de', 60 '#ffffff') # for log
set view map
set xrange [-10/Q0:10/Q0]
set yrange [-5/Q0:5/Q0]
set ylabel "{/Times-Italic Q}_c/{/Times-Italic Q}_0"
set xlabel "{/Times-Italic Q}_t/{/Times-Italic Q}_0"
#set cbrange [0.0:0.17] #nonlog
unset contour
set surface
set view map
splot "fig_8_excited_wavepacket.txt" u ($2/Q0):($1/Q0):($4 < 0.001 ? 1/0 : -log($4))  with pm3d at b notitle
set cntrparam levels incremental 0.120,0.005,0.20
set cntrlabel start 5 interval 100 onecolor
unset surface
set contour
set isosamples 500,500
splot f(x*Q0,y*Q0) with lines lc "black" notitle


unset multiplot


