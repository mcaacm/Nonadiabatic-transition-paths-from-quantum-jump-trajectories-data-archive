set terminal postscript enhanced colour font "Times-Roman,20"

# Conitcal intersection parameters
Q0 = (4.2538666725158691 + 4.8292193412780762)/2.0
kappa1 = -0.003858675*3.0
kappa2 = 0.0054756457*2.4
omegasc = 0.0043364174
omegast = 0.0027194482
e1 = 0.144792242
e2 = 0.177866612*0.87
lambda = 0.0096283166*0.3
a = 0.123646
# Lower and upper adiabatic surfaces of the conical intersection to plot heat maps
f(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 - sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))
 g(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 + sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2))

set output "figure_5.ps"
set multiplot

set label "|{/Symbol-Oblique f}_1><{/Symbol-Oblique f}_1|" at screen 0.27,0.63
set label "|{/Symbol-Oblique f}_2><{/Symbol-Oblique f}_2|" at screen 0.77,0.63
set label "a)" at screen 0.10,0.63 font "Times-Roman,24"
set label "b)" at screen 0.10,0.34 font "Times-Roman,24"

set ylabel "{/Times-Italic E}/eV"
set xlabel "Step"
set xrange [0:17]
set size 1.0,0.4
set origin 0.0,0.6
set yrange [0:0.65]
a = 0.123646*27.2113
set label 'a)' at 6.8,0.43 font "Times-Roman,24"
set label 'b)' at 7.8,0.43 font "Times-Roman,24"
plot "fig_5_pathway_energies.txt" u 1:($2*27.2113 - a) with lines lw 6 lc "black" notitle, "fig_5_pathway_energies.txt" u 1:($2*27.2113 -a) pt 5 ps 3 lc "black" notitle
unset label

unset colorbox
unset cbtics
unset xtics
unset xlabel

set size 0.55,0.4
set origin 0.05,0.25
set palette defined (20 '#053061', 30 '#2166ac', 40 '#4393c3', 50 '#92c5de', 60 '#ffffff') 
set view map
set xrange [-10/Q0:10/Q0]
set yrange [-5/Q0:5/Q0]
set ylabel "{/Times-Italic Q}_c/{/Times-Italic Q}_0"
set cbrange [0:3]
unset contour
set surface
set view map
splot "fig_5_pre_committor_wf.txt" u ($2/Q0):($1/Q0):($4 < 0.001 ? 1/0 : -log10($4))  with pm3d at b notitle
set cntrparam levels incremental 0.120,0.005,0.20
set cntrlabel start 5 interval 100 onecolor
unset surface
set contour
set isosamples 500,500
splot f(x*Q0,y*Q0) with lines lc "black" notitle

set origin 0.49,0.25
unset ylabel
unset cblabel
unset contour
unset ytics

set surface
set view map
#set palette defined (10 '#ffffff', 20 '#fddbc7', 30 '#f4a582', 40 '#d6604d', 50 '#b2182b', 60 '#67001f')  # For nonlog
set palette defined (20 '#67001f', 30 '#b2182b', 40 '#d6604d', 50 '#f4a582', 60 '#ffffff') #,110 '#053061') # For log
splot "fig_5_pre_committor_wf.txt" u ($2/Q0):($1/Q0):($3 < 0.001 ? 1/0 : -log10($3)) with pm3d at b notitle 
set cntrparam levels incremental 0.120,0.005,0.20
set cntrlabel start 5 interval 100 onecolor
unset surface
set contour
set isosamples 500,500
splot f(x*Q0,y*Q0) with lines lc "black" notitle


set xlabel "{/Times-Italic Q}_t/{/Times-Italic Q}_0"
set xtics
set ytics
set size 0.55,0.4
set origin 0.05,-0.02
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
splot "fig_5_post_committor_wf.txt" u ($2/Q0):($1/Q0):($4 < 0.001 ? 1/0 : -log10($4))  with pm3d at b notitle
set cntrparam levels incremental 0.120,0.005,0.20
set cntrlabel start 5 interval 100 onecolor
unset surface
set contour
set isosamples 500,500
splot f(x*Q0,y*Q0) with lines lc "black" notitle

set origin 0.49,-0.02
unset ylabel
unset cblabel
unset contour
unset ytics

set surface
set view map
#set palette defined (10 '#ffffff', 20 '#fddbc7', 30 '#f4a582', 40 '#d6604d', 50 '#b2182b', 60 '#67001f')  # For nonlog
set palette defined (20 '#67001f', 30 '#b2182b', 40 '#d6604d', 50 '#f4a582', 60 '#ffffff') #,110 '#053061') # For log
splot "fig_5_post_committor_wf.txt" u ($2/Q0):($1/Q0):($3 < 0.001 ? 1/0 : -log10($3)) with pm3d at b notitle 
set cntrparam levels incremental 0.120,0.005,0.20
set cntrlabel start 5 interval 100 onecolor
unset surface
set contour
set isosamples 500,500
splot f(x*Q0,y*Q0) with lines lc "black" notitle

unset multiplot

