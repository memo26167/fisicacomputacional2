#set terminal png

set terminal gif animate delay 4
#delay 4 means that there will be 0.01 × 4 seconds between frames in your animation 1/0.04 = 25 fps
set output "animacion.gif"

#set hidden3d


stats 'datos.dat' u 1:3 nooutput
tamaño=int(sqrt(STATS_records/STATS_blank))
set dgrid3d tamaño,tamaño qnorm 2

set xrange[STATS_min_x:STATS_max_x]
set yrange[STATS_min_x:STATS_max_x]
set zrange[STATS_min_y:STATS_max_y]


do for [i=1:STATS_blank-1]{
	start=(i-1)*(tamaño*tamaño)+1
	end=i*(tamaño*tamaño)
	splot "< awk 'NF!=0 { print $0 }' datos.dat" every ::start::end with lines
}