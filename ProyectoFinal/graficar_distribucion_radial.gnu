set terminal pngcairo enhanced size 1024, 1024  
set output 'distribucion_radial.png'

stats 'distribucion_radial_80.dat' u 1:2 nooutput
corte_g_80 = sprintf('{/Symbol r}{/Symbol s}^3 = 1.03  g({/Symbol s}) = %.3f ', STATS_max_y)

stats 'distribucion_radial_70.dat' u 1:2 nooutput
corte_g_70 = sprintf('{/Symbol r}{/Symbol s}^3 = 0.9    g({/Symbol s}) = %.3f ', STATS_max_y)

stats 'distribucion_radial_39.dat' u 1:2 nooutput
corte_g_39 = sprintf('{/Symbol r}{/Symbol s}^3 = 0.5    g({/Symbol s}) = %.3f ', STATS_max_y)

stats 'distribucion_radial_24.dat' u 1:2 nooutput
corte_g_24 = sprintf('{/Symbol r}{/Symbol s}^3 = 0.3    g({/Symbol s}) = %.3f ', STATS_max_y)

set style data lines

set title "Distribucion Radial"
set xlabel "r/{/Symbol s}"
set xrange [1:3]
set ylabel "N"
set grid

set origin 0,0
plot 'distribucion_radial_80.dat' using 1:2 ti corte_g_80, 'distribucion_radial_70.dat' using 1:2 ti corte_g_70, 'distribucion_radial_39.dat' using 1:2 ti corte_g_39, 'distribucion_radial_24.dat' using 1:2 ti corte_g_24


