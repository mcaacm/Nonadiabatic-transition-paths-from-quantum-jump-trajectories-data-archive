set terminal postscript portrait colour font "Times-Roman,18"
set output "figure_10.ps"

set label "a)" at screen 0.05,0.96
set label "b)" at screen 0.05,0.63

set multiplot
set xlabel "Step"
unset key
set size 1.0,0.35
set ytics 0.4,0.1,1.0
set origin 0,0.64
set ylabel "{/Times-Italic P_{R|L}}"
plot "fig_10_path_03.txt" u 1:17 with lines lw 4 lc rgb "#9e0142" title "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0=0.3", "fig_10_path_03.txt" u 1:17 ps 1 pt 5 lc rgb "#9e0142" notitle, "fig_10_path_08.txt" u 1:17 with lines lw 4 lc rgb "#fee08b" title "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0=0.8", "fig_10_path_08.txt" u 1:17 ps 1 pt 5 lc rgb "#fee08b" notitle, "fig_10_path_10.txt" u 1:17 with lines lw 4 lc rgb "#abdda4" title "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0=1.0", "fig_10_path_10.txt" u 1:17 ps 1 pt 5 lc rgb "#abdda4" notitle, "fig_10_path_13.txt" u 1:17  with lines lw 4 lc rgb "#3288bd" title "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0=1.3", "fig_10_path_13.txt" u 1:17 ps 1 pt 5 lc rgb "#3288bd" notitle
set key horizontal at screen 0.85,0.25
set size 1.0,0.36
set origin 0,0.26
set ylabel "{/Times-Italic P_{L|R}}"
plot "fig_10_path_03_B.txt" u 1:17  with lines lw 4 lc rgb "#9e0142" title "{/Symbol-Oblique l}=0.3{/Symbol-Oblique l}_0", "fig_10_path_03_B.txt" u 1:17 ps 1 pt 5 lc rgb "#9e0142" notitle, "fig_10_path_08_B.txt" u 1:17  with lines lw 4 lc rgb "#fee08b" title "{/Symbol-Oblique l}=0.8{/Symbol-Oblique l}_0", "fig_10_path_08_B.txt" u 1:17 ps 1 pt 5 lc rgb "#fee08b" notitle, "fig_10_path_10_B.txt" u 1:17  with lines lw 4 lc rgb "#abdda4" title "{/Symbol-Oblique l}=1.0{/Symbol-Oblique l}_0", "fig_10_path_10_B.txt" u 1:17 ps 1 pt 5 lc rgb "#abdda4" notitle, "fig_10_path_13_B.txt" u 1:17  with lines lw 4 lc rgb "#3288bd" title "{/Symbol-Oblique l}=1.3{/Symbol-Oblique l}_0", "fig_10_path_13_B.txt" u 1:17 ps 1 pt 5 lc rgb "#3288bd" notitle

unset multiplot

