set terminal postscript enhanced colour font "Times-Roman,20"

set palette defined (20 '#67001f', 30 '#b2182b', 40 '#d6604d', 50 '#f4a582', 60 '#fddbc7', 70 '#d1e5f0', 80 '#92c5de', 90 '#4393c3', 100 '#2166ac',110 '#053061')

# System parameters
Q0 = (4.1687483787536621 + 4.7293429374694824) / 2.0  # Defining the relative coordinate Q_0
a = 0.12285201251506805*27.211  # Defining the zero of energy as that of the system's lowest energy eigenvector
# Conical intersection specification
kappa1 = -0.003858675*3.0 
kappa2 = 0.0054756457*2.4
omegasc = 0.0043364174
omegast = 0.0027194482
e1 = 0.144792242
e2 = 0.177866612*0.87
lambda = 0.0096283166*1.3

# The functions representing the two adiabatic potential energy surfaces as a function of Q_t (x) and Q_c (y)
f(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 - sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2)) 
g(x,y) = 0.5*(1.*e1 + 1.*e2 + 1.*kappa1*x + 1.*kappa2*x + 1.*omegast*x**2 + 1.*omegasc*y**2 + sqrt(1.*e1**2 - 2.*e1*e2 + 1.*e2**2 + 2.*e1*kappa1*x - 2.*e2*kappa1*x - 2.*e1*kappa2*x + 2.*e2*kappa2*x + 1.*kappa1**2*x**2 - 2.*kappa1*kappa2*x**2 + 1.*kappa2**2*x**2 + 4.*lambda**2*y**2)) 


# Q_t coordinate of the intersection. The intersection Q_c coordinate is 0.
intersect2 = (0.5*(-2.*e1*kappa1 + 2.*e2*kappa1 + 2.*e1*kappa2 - 2.*e2*kappa2 + sqrt((2.*e1*kappa1 - 2.*e2*kappa1 - 2.*e1*kappa2 + 2.*e2*kappa2)**2 - 4.*(1.*e1**2 - 2.*e1*e2 + 1.*e2**2)*(kappa1**2 - 2.*kappa1*kappa2 + kappa2**2))))/(kappa1**2 - 2.*kappa1*kappa2 + kappa2**2)


set output "figure_1.ps"
set multiplot

set label "a)" at screen 0.05,0.9
set label "b)" at screen 0.08,0.40
set label "c)" at screen 0.48,0.40

set size 0.9,0.8
set origin 0.05,0.4

unset xrange
unset yrange
unset zrange
unset cbrange
# Label axes
set xlabel "Q_{t}/Q_0" font "Times-Italic,20"
set ylabel "Q_{c}/Q_0" font "Times-Italic,20"
set zrange[3.2 - a:5.5 - a]
set cbrange[3.2 - a:5.5 - a]
# Set location of 3d plat base plane
set xyplane at 3.2-a
# Adjust axes to 0
set yrange [-4/Q0:4/Q0]
set xrange [-7/Q0:7/Q0]
set hidden3d
set ytics -0.8,0.4,0.8

set zlabel "{/Times-Italic E}/eV" rotate parallel

set view 50,30,1,1  # Rotate plane of view
set isosamples 800  # Sample enough points for smooth plotting
set samples 1000
# Plot the adiabatic potential energy surfaces. Plot energies relative to the zero of energy at 'a'.
splot 27.211*f(x*Q0,y*Q0) - a with pm3d notitle, 27.211*g(x*Q0,y*Q0) - a with pm3d notitle 

set size 0.45,0.6
set origin 0.05,-0.15
set yrange [-0.1:1.0]
set ytics 0,0.3,1.2
set xtics -1.5,0.5,1.5
set xrange [-7/Q0:7/Q0]
set ylabel "{/Times-Italic E}/eV" font "Times-Roman,20"
set xlabel "Q_{t}/Q_0" 
# Plot the locations of eigenvectors and a cut through the lower adiabatic surface at the conical intersection (where Q_c is 0)
plot "fig_1_eigenvectors.txt" u ($2/Q0):(($3)*27.211 - a) with lines lw 1 lc "grey" notitle, 27.211*(f(x*Q0,0)) - a with lines lw 3 lc rgb "#9e0142" notitle, 27.211*(g(x*Q0,0)) - a with lines lw 3 lc rgb "#3288bd" notitle 
set xlabel "Q_{c}/Q_0" 
#set xtics -4/Q0,4/Q0,4/Q0
set xtics -0.8,0.4,0.8
set size 0.45,0.6
set origin 0.47,-0.15
set xrange [-4/Q0:4/Q0]
# Plot the locations of eigenvectors and a cut through the lower adiabatic surface at the conical intersection (where Q_t=instersect2)
plot "fig_1_eigenvectors_4.txt" u ($2/Q0):(($3)*27.211 - a) with lines lw 1 lc "grey" notitle, 27.211*(f(intersect2,x*Q0)) - a with lines lw 3 lc rgb "#9e0142" notitle, 27.211*(g(intersect2,x*Q0)) - a with lines lw 3 lc rgb "#3288bd" notitle 

unset multiplot
