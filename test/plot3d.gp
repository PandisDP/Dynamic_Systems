set contour base 
set xrange [0:200] 
set term postscript color  10 
set output "contour.ps"
set xlabel 'Time(Units)' 
set ylabel 'Outcome Value [Min-Max]' 
set zlabel 'Count'
set hidden3d 
set cntrparam levels 8
set ticslevel 1
set title "Frequency values of: Poblacion"
splot 'contour.0001' t '' w l
set title "Frequency values of: Minas"
splot 'contour.0002' t '' w l
set title "Frequency values of: Ingresos"
splot 'contour.0003' t '' w l
set title "Frequency values of: Deudas"
splot 'contour.0004' t '' w l
set title "Frequency values of: Presupuesto"
splot 'contour.0005' t '' w l
set title "Frequency values of: PBIpc"
splot 'contour.0006' t '' w l
set title "Frequency values of: DEUpc"
splot 'contour.0007' t '' w l
exit
