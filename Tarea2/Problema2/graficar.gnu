set terminal pngcairo enhanced size 1024, 1024 
set output 'grafico_solucion.png'
set style data lines 
set title "Solucion"
set xlabel "Eje x"
set ylabel "Eje y"
set grid
plot 'datos.dat' using 1:2 ti "Solucion"
