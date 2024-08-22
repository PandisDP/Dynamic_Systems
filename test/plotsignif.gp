set term postscript color  10 
set output "signif.ps"
set yrange [0:0.3] 
set ytics 0.01 
set xrange [0:350] 
set xlabel 'Time(Units)' 
set title "PRCC Significance: trucks_in"
plot 'prcc.0001' u 1:10 t 'trucks_init', 'prcc.0001' u 1:12 t 'percentage_dis', 'prcc.0001' u 1:14 t 'distance_speed', 'prcc.0001' u 1:16 t 'flow_left', 'prcc.0001' u 1:18 t 'flow_right', 'prcc.0001' u 1:20 t 'flow_up', 'prcc.0001' u 1:22 t 'flow_down', 'prcc.0001' u 1:24 t 'percentage_minimum_dist', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: queue_line_a"
plot 'prcc.0002' u 1:10 t 'trucks_init', 'prcc.0002' u 1:12 t 'percentage_dis', 'prcc.0002' u 1:14 t 'distance_speed', 'prcc.0002' u 1:16 t 'flow_left', 'prcc.0002' u 1:18 t 'flow_right', 'prcc.0002' u 1:20 t 'flow_up', 'prcc.0002' u 1:22 t 'flow_down', 'prcc.0002' u 1:24 t 'percentage_minimum_dist', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: queue_line_b"
plot 'prcc.0003' u 1:10 t 'trucks_init', 'prcc.0003' u 1:12 t 'percentage_dis', 'prcc.0003' u 1:14 t 'distance_speed', 'prcc.0003' u 1:16 t 'flow_left', 'prcc.0003' u 1:18 t 'flow_right', 'prcc.0003' u 1:20 t 'flow_up', 'prcc.0003' u 1:22 t 'flow_down', 'prcc.0003' u 1:24 t 'percentage_minimum_dist', 0 t '' w l, 0.01 t '' w l 
exit
