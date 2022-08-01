set terminal postscript enhanced colour font "Times-Roman,18"

set output "figure_4.ps"
#set label "a)" at screen 0.05,0.96
#set label "b)" at screen 0.05,0.63
#set label "c)" at screen 0.05,0.30

set key bottom right
set xlabel "Step"
set ylabel "{/Times-Italic P_{L|R}}"
set yrange [-0.01:1.01]
set size 1.0,0.8
plot "fig_4_path_03.txt" u 1:17 with lines lw 4 lc rgb "#9e0142" title "{/Symbol-Oblique l}=0.3{/Symbol-Oblique l}_0", "fig_4_path_03.txt" u 1:17 ps 1 pt 5 lc rgb "#9e0142" notitle, "fig_4_path_08.txt" u 1:17 with lines lw 4 lc rgb "#fee08b" title "{/Symbol-Oblique l}=0.8{/Symbol-Oblique l}_0", "fig_4_path_08.txt" u 1:17 ps 1 pt 5 lc rgb "#fee08b" notitle, "fig_4_path_10.txt" u 1:17 with lines lw 4 lc rgb "#abdda4" title "{/Symbol-Oblique l}=1.0{/Symbol-Oblique l}_0", "fig_4_path_10.txt" u 1:17 ps 1 pt 5 lc rgb "#abdda4" notitle, "fig_4_path_13.txt" u 1:17 with lines lw 4 lc rgb "#3288bd" title "{/Symbol-Oblique l}=1.3{/Symbol-Oblique l}_0", "fig_4_path_13.txt" u 1:17 ps 1 pt 5 lc rgb "#3288bd" notitle, "fig_4_committor_points.txt" u 1:17 ps 1 pt 35 lc "black" notitle



