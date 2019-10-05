set terminal pngcairo enhanced size 800, 800 
set output 'grafico_solucion.png'
set style data lines 
set title "Solucion"
set xlabel "Eje x"
set ylabel "Eje y"
set grid
plot for [IDX=0:4] 'salida.dat' i IDX u 1:2 w lines title columnheader(1)
