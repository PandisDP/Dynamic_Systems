set term postscript color  10 
set output "stats.ps"
set xlabel 'Time(Units)' 
set title "Stats :Poblacion"
plot 'prcc.0001' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0001' u 1:6 t 'Min' w l, 'prcc.0001' u 1:5 t 'Max' w l
set title "Stats :Minas"
plot 'prcc.0002' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0002' u 1:6 t 'Min' w l, 'prcc.0002' u 1:5 t 'Max' w l
set title "Stats :Ingresos"
plot 'prcc.0003' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0003' u 1:6 t 'Min' w l, 'prcc.0003' u 1:5 t 'Max' w l
set title "Stats :Deudas"
plot 'prcc.0004' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0004' u 1:6 t 'Min' w l, 'prcc.0004' u 1:5 t 'Max' w l
set title "Stats :Presupuesto"
plot 'prcc.0005' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0005' u 1:6 t 'Min' w l, 'prcc.0005' u 1:5 t 'Max' w l
set title "Stats :PBIpc"
plot 'prcc.0006' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0006' u 1:6 t 'Min' w l, 'prcc.0006' u 1:5 t 'Max' w l
set title "Stats :DEUpc"
plot 'prcc.0007' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0007' u 1:6 t 'Min' w l, 'prcc.0007' u 1:5 t 'Max' w l
exit
