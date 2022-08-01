set terminal postscript colour font "Times-Roman,15"
set output "figure_2.ps"

set multiplot

set label 'a)' at screen 0.18,0.67
set label 'b)' at screen 0.52,0.67

set origin 0.15,0.2
set size 0.35,0.5  

set font "Times-Roman,24"
set xrange [1:5000]
set key top left
set xlabel "{/Times-Italic t}/fs"
set ylabel "{/Symbol-Oblique r}_{/Times-Italic R}"
set logscale 
# The conversion factor of 0.024188 is required for atomic unit to fs 
plot "fig_2_density_matrix_03.txt" u ($2*0.024188):3 with lines lw 4 lc rgb "#9e0142" title "{/Symbol-Oblique l}=0.3{/Symbol-Oblique l}_0", "fig_2_density_matrix_13.txt" u ($2*0.024188):3 with lines lw 4 lc rgb "#3288bd" title "{/Symbol-Oblique l}=1.3{/Symbol-Oblique l}_0"
unset logscale

set origin 0.48,0.17
set size 0.35,0.53  
set xlabel "{/Symbol-Oblique l}/{/Symbol-Oblique l}_0"
unset xrange
set logscale y
set xtics 0.0,0.4,1.6
set ylabel "{/Times-Italic k_{L,R}}/fs"
plot "fig_2_flux_Bto1.txt" u 1:($2/0.0241888) with lines lw 4 lc "black" notitle, "fig_2_flux_Bto1.txt" u 1:($2/0.0241888) ps 1 pt 5 lc "black" notitle
unset logscale

unset multiplot


