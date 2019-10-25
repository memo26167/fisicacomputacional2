set terminal pngcairo enhanced size 1024, 1024 
set output 'grafico_solucion.png'
set style data lines 
set title "Solucion"
set xlabel "Eje x"
set yrange [0:1.5]
set ylabel "Eje y"
set grid
set style line 1 lt rgb "red" lw 3
set style line 2 lt rgb "yellow" lw 2
set style line 3 lt rgb "blue" lw 3

plot 'datos.dat' using 1:2 ls 1 ti "FTBS" ,'datos.dat' using 1:3 ls 2 ti "FTCS",'datos.dat' using 1:4 ls 3 ti "LAX"
