set term postscript color  10 
set output "stats.ps"
set xlabel 'Time(Units)' 
set title "Stats :BE"
plot 'prcc.0001' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0001' u 1:6 t 'Min' w l, 'prcc.0001' u 1:5 t 'Max' w l
set title "Stats :T0"
plot 'prcc.0002' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0002' u 1:6 t 'Min' w l, 'prcc.0002' u 1:5 t 'Max' w l
set title "Stats :T1"
plot 'prcc.0003' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0003' u 1:6 t 'Min' w l, 'prcc.0003' u 1:5 t 'Max' w l
set title "Stats :T2"
plot 'prcc.0004' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0004' u 1:6 t 'Min' w l, 'prcc.0004' u 1:5 t 'Max' w l
set title "Stats :T8"
plot 'prcc.0005' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0005' u 1:6 t 'Min' w l, 'prcc.0005' u 1:5 t 'Max' w l
set title "Stats :MR"
plot 'prcc.0006' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0006' u 1:6 t 'Min' w l, 'prcc.0006' u 1:5 t 'Max' w l
set title "Stats :MI"
plot 'prcc.0007' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0007' u 1:6 t 'Min' w l, 'prcc.0007' u 1:5 t 'Max' w l
set title "Stats :MA"
plot 'prcc.0008' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0008' u 1:6 t 'Min' w l, 'prcc.0008' u 1:5 t 'Max' w l
set title "Stats :IG"
plot 'prcc.0009' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0009' u 1:6 t 'Min' w l, 'prcc.0009' u 1:5 t 'Max' w l
set title "Stats :I12"
plot 'prcc.0010' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0010' u 1:6 t 'Min' w l, 'prcc.0010' u 1:5 t 'Max' w l
set title "Stats :I4"
plot 'prcc.0011' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0011' u 1:6 t 'Min' w l, 'prcc.0011' u 1:5 t 'Max' w l
set title "Stats :I10"
plot 'prcc.0012' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0012' u 1:6 t 'Min' w l, 'prcc.0012' u 1:5 t 'Max' w l
set title "Stats :FA"
plot 'prcc.0013' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0013' u 1:6 t 'Min' w l, 'prcc.0013' u 1:5 t 'Max' w l
set title "Stats :BI"
plot 'prcc.0014' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0014' u 1:6 t 'Min' w l, 'prcc.0014' u 1:5 t 'Max' w l
set title "Stats :TB"
plot 'prcc.0015' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0015' u 1:6 t 'Min' w l, 'prcc.0015' u 1:5 t 'Max' w l
set title "Stats :T"
plot 'prcc.0016' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0016' u 1:6 t 'Min' w l, 'prcc.0016' u 1:5 t 'Max' w l
set title "Stats :M"
plot 'prcc.0017' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0017' u 1:6 t 'Min' w l, 'prcc.0017' u 1:5 t 'Max' w l
set title "Stats :MIp"
plot 'prcc.0018' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0018' u 1:6 t 'Min' w l, 'prcc.0018' u 1:5 t 'Max' w l
set title "Stats :MITIratio"
plot 'prcc.0019' u 1:2:7:8 t 'Avg' w errorbars, 'prcc.0019' u 1:6 t 'Min' w l, 'prcc.0019' u 1:5 t 'Max' w l
exit
