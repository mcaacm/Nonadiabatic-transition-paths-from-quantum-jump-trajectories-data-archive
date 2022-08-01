set terminal postscript enhanced colour font "Times-Roman,25"

set palette defined (0 "green", 1 "blue", 2 "red", 3 "orange" )

set output "figure_3.ps"
set key bottom right
set xrange [0:40]
set size 1.0,0.8
set xlabel "{/Times-Italic n}" font "Times-Roman,25"
set ylabel "{/Times-Italic F_n}" font "Times-Roman,25"
plot "fig_3_paths_03.txt" u 1:3 with lines lw 7 lc rgb "#9e0142" title "{/Symbol-Oblique l}=0.3{/Symbol-Oblique l}_0", "fig_3_paths_05.txt" u 1:3 with lines lw 7 lc rgb "#f46d43" title "{/Symbol-Oblique l}=0.5{/Symbol-Oblique l}_0", "fig_3_paths_08.txt" u 1:3 with lines lw 7 lc rgb "#fee08b" title "{/Symbol-Oblique l}=0.8{/Symbol-Oblique l}_0", "fig_3_paths_10.txt" u 1:3 with lines lw 7 lc rgb "#abdda4" title "{/Symbol-Oblique l}=1.0{/Symbol-Oblique l}_0", "fig_3_paths_13.txt" u 1:3 with lines lw 7 lc rgb "#3288bd"  title "{/Symbol-Oblique l}=1.3{/Symbol-Oblique l}_0", "fig_3_paths_15.txt" u 1:3 with lines lw 7 lc rgb "#5e4fa2" title "{/Symbol-Oblique l}=1.5{/Symbol-Oblique l}_0"


