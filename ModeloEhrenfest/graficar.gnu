set terminal pngcairo enhanced size 2048, 1024 
set output 'evolucion.png'
set style data lines 
set title "Evolucion"
set xlabel "Tiempo"
set ylabel "Numero de moleculas"
set grid
plot 'evolucion.dat' every :::0::0 ti "Evolucion" with steps

set output 'distribucion.png'
set style data lines 
set title "Distribucion"
set xlabel "Numero de moleculas"
set ylabel "Probabilidad"
set grid
plot 'evolucion.dat' every :::1::1 ti "Distribucion"