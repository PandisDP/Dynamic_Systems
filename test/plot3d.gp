set contour base 
set xrange [0:225] 
set term postscript color  10 
set output "contour.ps"
set xlabel 'Time(Units)' 
set ylabel 'Outcome Value [Min-Max]' 
set zlabel 'Count'
set hidden3d 
set cntrparam levels 8
set ticslevel 1
set title "Frequency values of: variable1"
splot 'contour.0001' t '' w l
set title "Frequency values of: variable2"
splot 'contour.0002' t '' w l
exit
