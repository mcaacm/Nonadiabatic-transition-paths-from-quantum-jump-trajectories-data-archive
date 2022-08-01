set terminal postscript enhanced colour font "Times-Roman,20"
set output "figure_9.ps"
set size 1.0,0.8
set xrange [0:200]
set ylabel "{/Times-Italic F_n}"
set key bottom right
set xlabel "{/Times-Italic n}"
plot "fig_9_cflux_03.txt" u 1:3 with lines lw 7 lc rgb "#9e0142" title "{/Symbol-Oblique l}=0.3{/Symbol-Oblique l}_0", "fig_9_cflux_05.txt" u 1:3 with lines lw 7 lc "#f46d43" title "{/Symbol-Oblique l}=0.5{/Symbol-Oblique l}_0","fig_9_cflux_08.txt" u 1:3 with lines lw 7 lc rgb "#fee08b" title "{/Symbol-Oblique l}=0.8{/Symbol-Oblique l}_0", "fig_9_cflux_10.txt" u 1:3 with lines lw 7 lc rgb "#abdda4" title "{/Symbol-Oblique l}=1.0{/Symbol-Oblique l}_0", "fig_9_cflux_13.txt" u 1:3 with lines lw 7 lc rgb "#3288bd" title "{/Symbol-Oblique l}=1.3{/Symbol-Oblique l}_0", "fig_9_cflux_15.txt" u 1:3 with lines lw 7 lc "#5e4fa2" title "{/Symbol-Oblique l}=1.5{/Symbol-Oblique l}_0"



