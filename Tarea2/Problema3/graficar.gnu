set terminal pngcairo enhanced size 1024, 1024 
set output 'grafico_solucion.png'
set style data lines 
set title "Solucion"
set xlabel "Eje x"
set yrange [0:1.5]
set ylabel "Eje y"
set grid
set style line 1 lt rgb "red" lw 3
set style line 2 lt 1 lc 2 lw 5
set style line 3 lt 2 lc 3 lw 2
set style line 4 lt rgb "black" lw 2

plot 'datos.dat' using 1:2 ls 1 ti "FTBS" ,'datos.dat' using 1:3 ls 2 ti "FTCS", 'datos.dat' using 1:4 ls 3 ti "BTCS",'datos.dat' using 1:5 ls 4 ti "LAX"
