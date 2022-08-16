set terminal postscript enhanced colour landscape font "Times-Roman,28"
set output "figure_13.ps"

set ytics 0.1
set ylabel "({/Times-Italic E - E_0})/eV"
set xlabel "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0"
plot "fig_13_postc_avg_vertical.txt" u 1:(($2-$8)*27.211) with lines lw 5 lc rgb "#3288bd" notitle, "fig_13_postc_avg_vertical.txt" u 1:(($2-$8)*27.211) ps 2 pt 5 lc rgb "#3288bd" notitle, "fig_13_coupling_hte.txt" u 1:(($2-$5)*27.211) with lines lw 5 dt "-" lc "black" title "Barrier Energy"


