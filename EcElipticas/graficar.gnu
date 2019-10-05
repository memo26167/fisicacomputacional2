set terminal pngcairo enhanced size 1024, 1024 
set output 'grafico_solucion.png'
set title "Solucion"
set pm3d
set xlabel "Eje x"
set ylabel "Eje y"
set zlabel "Eje z"
set grid
splot 'salida.dat' with pm3d ti "Solucion"
