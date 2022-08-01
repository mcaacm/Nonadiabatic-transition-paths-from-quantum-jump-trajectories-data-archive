set terminal postscript enhanced colour font "Times-Roman,24"

set key bottom left
set output "figure_7.ps"
set xlabel "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0"
set ylabel "{/Times-Italic (E - E_0)}/eV"
plot "fig_7_postc_averages.txt" u 1:(($2-$8)*27.211) with lines lw 6 lc rgb "#3288bd" title "Post-Jump", "fig_7_prec_averages.txt" u 1:(($2-$8)*27.211) with lines lw 6 lc rgb "#9e0142" title "Pre-Jump", "fig_7_coupling_hte.txt" u 1:(($2-$5)*27.2113) with lines lw 6 dt "." lc "black" title "Barrier Energy"

