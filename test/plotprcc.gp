set term postscript color  10 
set output "prccvals.ps"
set ytics 0.2 
set yrange [-1:1] 
set xrange [0:500] 
set xlabel 'Time(Units)' 
set title "PRCC : BE"
plot 'prcc.0001' u 1:9 t 'alpha2', 'prcc.0001' u 1:11 t 'alpha4', 'prcc.0001' u 1:13 t 'alpha3', 'prcc.0001' u 1:15 t 'alpha10', 'prcc.0001' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : T0"
plot 'prcc.0002' u 1:9 t 'alpha2', 'prcc.0002' u 1:11 t 'alpha4', 'prcc.0002' u 1:13 t 'alpha3', 'prcc.0002' u 1:15 t 'alpha10', 'prcc.0002' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : T1"
plot 'prcc.0003' u 1:9 t 'alpha2', 'prcc.0003' u 1:11 t 'alpha4', 'prcc.0003' u 1:13 t 'alpha3', 'prcc.0003' u 1:15 t 'alpha10', 'prcc.0003' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : T2"
plot 'prcc.0004' u 1:9 t 'alpha2', 'prcc.0004' u 1:11 t 'alpha4', 'prcc.0004' u 1:13 t 'alpha3', 'prcc.0004' u 1:15 t 'alpha10', 'prcc.0004' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : T8"
plot 'prcc.0005' u 1:9 t 'alpha2', 'prcc.0005' u 1:11 t 'alpha4', 'prcc.0005' u 1:13 t 'alpha3', 'prcc.0005' u 1:15 t 'alpha10', 'prcc.0005' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : MR"
plot 'prcc.0006' u 1:9 t 'alpha2', 'prcc.0006' u 1:11 t 'alpha4', 'prcc.0006' u 1:13 t 'alpha3', 'prcc.0006' u 1:15 t 'alpha10', 'prcc.0006' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : MI"
plot 'prcc.0007' u 1:9 t 'alpha2', 'prcc.0007' u 1:11 t 'alpha4', 'prcc.0007' u 1:13 t 'alpha3', 'prcc.0007' u 1:15 t 'alpha10', 'prcc.0007' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : MA"
plot 'prcc.0008' u 1:9 t 'alpha2', 'prcc.0008' u 1:11 t 'alpha4', 'prcc.0008' u 1:13 t 'alpha3', 'prcc.0008' u 1:15 t 'alpha10', 'prcc.0008' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : IG"
plot 'prcc.0009' u 1:9 t 'alpha2', 'prcc.0009' u 1:11 t 'alpha4', 'prcc.0009' u 1:13 t 'alpha3', 'prcc.0009' u 1:15 t 'alpha10', 'prcc.0009' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : I12"
plot 'prcc.0010' u 1:9 t 'alpha2', 'prcc.0010' u 1:11 t 'alpha4', 'prcc.0010' u 1:13 t 'alpha3', 'prcc.0010' u 1:15 t 'alpha10', 'prcc.0010' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : I4"
plot 'prcc.0011' u 1:9 t 'alpha2', 'prcc.0011' u 1:11 t 'alpha4', 'prcc.0011' u 1:13 t 'alpha3', 'prcc.0011' u 1:15 t 'alpha10', 'prcc.0011' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : I10"
plot 'prcc.0012' u 1:9 t 'alpha2', 'prcc.0012' u 1:11 t 'alpha4', 'prcc.0012' u 1:13 t 'alpha3', 'prcc.0012' u 1:15 t 'alpha10', 'prcc.0012' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : FA"
plot 'prcc.0013' u 1:9 t 'alpha2', 'prcc.0013' u 1:11 t 'alpha4', 'prcc.0013' u 1:13 t 'alpha3', 'prcc.0013' u 1:15 t 'alpha10', 'prcc.0013' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : BI"
plot 'prcc.0014' u 1:9 t 'alpha2', 'prcc.0014' u 1:11 t 'alpha4', 'prcc.0014' u 1:13 t 'alpha3', 'prcc.0014' u 1:15 t 'alpha10', 'prcc.0014' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : TB"
plot 'prcc.0015' u 1:9 t 'alpha2', 'prcc.0015' u 1:11 t 'alpha4', 'prcc.0015' u 1:13 t 'alpha3', 'prcc.0015' u 1:15 t 'alpha10', 'prcc.0015' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : T"
plot 'prcc.0016' u 1:9 t 'alpha2', 'prcc.0016' u 1:11 t 'alpha4', 'prcc.0016' u 1:13 t 'alpha3', 'prcc.0016' u 1:15 t 'alpha10', 'prcc.0016' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : M"
plot 'prcc.0017' u 1:9 t 'alpha2', 'prcc.0017' u 1:11 t 'alpha4', 'prcc.0017' u 1:13 t 'alpha3', 'prcc.0017' u 1:15 t 'alpha10', 'prcc.0017' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : MIp"
plot 'prcc.0018' u 1:9 t 'alpha2', 'prcc.0018' u 1:11 t 'alpha4', 'prcc.0018' u 1:13 t 'alpha3', 'prcc.0018' u 1:15 t 'alpha10', 'prcc.0018' u 1:17 t 'k18', 0 t '' w l 
set title "PRCC : MITIratio"
plot 'prcc.0019' u 1:9 t 'alpha2', 'prcc.0019' u 1:11 t 'alpha4', 'prcc.0019' u 1:13 t 'alpha3', 'prcc.0019' u 1:15 t 'alpha10', 'prcc.0019' u 1:17 t 'k18', 0 t '' w l 
exit
