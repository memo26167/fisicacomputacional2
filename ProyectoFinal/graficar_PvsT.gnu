set terminal pngcairo enhanced size 1024, 1024  
set output 'PresionVsTemperatura.png'

pi = 3.14159265359

stats 'distribucion_radial_80.dat' u 1:2 nooutput
corte_g_80 = STATS_max_y
m_80 = (1 + 4*pi* 80 *0.012788868 / 6 * corte_g_80) * 80 * 1.38064852e-23

stats 'distribucion_radial_70.dat' u 1:2 nooutput
corte_g_70 = STATS_max_y
m_70 = (1 + 4*pi* 70 *0.012788868 / 6 * corte_g_70) * 70 * 1.38064852e-23

stats 'distribucion_radial_39.dat' u 1:2 nooutput
corte_g_39 = STATS_max_y
m_39 = (1 + 4*pi* 39 *0.012788868 / 6 * corte_g_39) * 39 * 1.38064852e-23

stats 'distribucion_radial_24.dat' u 1:2 nooutput
corte_g_24 = STATS_max_y
m_24 = (1 + 4*pi* 24 *0.012788868 / 6 * corte_g_24) * 24 * 1.38064852e-23

set style data lines

set title "Presion v/s Temperatura"
set xlabel "Temperatura"
set xrange [0:300]
set ylabel "Presion"
set grid

plot m_80*x ti '{/Symbol r}{/Symbol s}^3 = 1.03', m_70*x ti '{/Symbol r}{/Symbol s}^3 = 0.9', m_39*x ti '{/Symbol r}{/Symbol s}^3 = 0.5', m_24*x ti '{/Symbol r}{/Symbol s}^3 = 0.3'

