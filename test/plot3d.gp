set contour base 
set xrange [0:500] 
set term postscript color  10 
set output "contour.ps"
set xlabel 'Time(Units)' 
set ylabel 'Outcome Value [Min-Max]' 
set zlabel 'Count'
set hidden3d 
set cntrparam levels 8
set ticslevel 1
set title "Frequency values of: BE"
splot 'contour.0001' t '' w l
set title "Frequency values of: T0"
splot 'contour.0002' t '' w l
set title "Frequency values of: T1"
splot 'contour.0003' t '' w l
set title "Frequency values of: T2"
splot 'contour.0004' t '' w l
set title "Frequency values of: T8"
splot 'contour.0005' t '' w l
set title "Frequency values of: MR"
splot 'contour.0006' t '' w l
set title "Frequency values of: MI"
splot 'contour.0007' t '' w l
set title "Frequency values of: MA"
splot 'contour.0008' t '' w l
set title "Frequency values of: IG"
splot 'contour.0009' t '' w l
set title "Frequency values of: I12"
splot 'contour.0010' t '' w l
set title "Frequency values of: I4"
splot 'contour.0011' t '' w l
set title "Frequency values of: I10"
splot 'contour.0012' t '' w l
set title "Frequency values of: FA"
splot 'contour.0013' t '' w l
set title "Frequency values of: BI"
splot 'contour.0014' t '' w l
set title "Frequency values of: TB"
splot 'contour.0015' t '' w l
set title "Frequency values of: T"
splot 'contour.0016' t '' w l
set title "Frequency values of: M"
splot 'contour.0017' t '' w l
set title "Frequency values of: MIp"
splot 'contour.0018' t '' w l
set title "Frequency values of: MITIratio"
splot 'contour.0019' t '' w l
exit
