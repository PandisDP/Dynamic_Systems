set contour base 
set xrange [0:3600] 
set term postscript color  10 
set output "contour.ps"
set xlabel 'Time(Units)' 
set ylabel 'Outcome Value [Min-Max]' 
set zlabel 'Count'
set hidden3d 
set cntrparam levels 8
set ticslevel 1
set title "Frequency values of: trucks_in"
splot 'contour.0001' t '' w l
set title "Frequency values of: queue_line_a"
splot 'contour.0002' t '' w l
set title "Frequency values of: queue_line_b"
splot 'contour.0003' t '' w l
set title "Frequency values of: flow_trucks"
splot 'contour.0004' t '' w l
exit
