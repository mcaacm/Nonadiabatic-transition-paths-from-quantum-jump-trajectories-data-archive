set terminal postscript colour font "Times-Roman,22"

set output "figure_15.ps"
const = 10**-16
a = 0.02418884254
set xlabel "{/Symbol-Oblique t}/fs"
set ylabel "{/Times-Italic k_{R,L}}10^{-16}/fs"
set size 1.0,0.8
set logscale x
plot "fig_15_flux_vs_tau.txt" u ($2*a):($7/(a*const)) with lines lw 7 lc "black" notitle, "fig_15_flux_vs_tau.txt" u ($2*a):($7/(a*const)) pt 5 ps 2 lc "black" notitle


