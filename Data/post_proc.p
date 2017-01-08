set autoscale
set xr [0:6.28]
set xlabel " x "
set ylabel " u(x,t) "

plot "physical_init.dat"  using 1:2 with lines title "t = 0s", \
	 "physical_final.dat" using 1:2 with lines title "t_final"

