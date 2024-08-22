set term postscript color  10 
set output "prccvals.ps"
set ytics 0.2 
set yrange [-1:1] 
set xrange [0:3600] 
set xlabel 'Time(Units)' 
set title "PRCC : trucks_in"
plot 'prcc.0001' u 1:9 t 'trucks_init', 'prcc.0001' u 1:11 t 'percentage_dis', 'prcc.0001' u 1:13 t 'distance_speed', 'prcc.0001' u 1:15 t 'flow_left', 'prcc.0001' u 1:17 t 'flow_right', 'prcc.0001' u 1:19 t 'flow_up', 'prcc.0001' u 1:21 t 'flow_down', 'prcc.0001' u 1:23 t 'percentage_minimum_dist', 0 t '' w l 
set title "PRCC : queue_line_a"
plot 'prcc.0002' u 1:9 t 'trucks_init', 'prcc.0002' u 1:11 t 'percentage_dis', 'prcc.0002' u 1:13 t 'distance_speed', 'prcc.0002' u 1:15 t 'flow_left', 'prcc.0002' u 1:17 t 'flow_right', 'prcc.0002' u 1:19 t 'flow_up', 'prcc.0002' u 1:21 t 'flow_down', 'prcc.0002' u 1:23 t 'percentage_minimum_dist', 0 t '' w l 
set title "PRCC : queue_line_b"
plot 'prcc.0003' u 1:9 t 'trucks_init', 'prcc.0003' u 1:11 t 'percentage_dis', 'prcc.0003' u 1:13 t 'distance_speed', 'prcc.0003' u 1:15 t 'flow_left', 'prcc.0003' u 1:17 t 'flow_right', 'prcc.0003' u 1:19 t 'flow_up', 'prcc.0003' u 1:21 t 'flow_down', 'prcc.0003' u 1:23 t 'percentage_minimum_dist', 0 t '' w l 
exit
