set term postscript color  10 
set output "signif.ps"
set yrange [0:0.3] 
set ytics 0.01 
set xrange [0:500] 
set xlabel 'Time(Units)' 
set title "PRCC Significance: BE"
plot 'prcc.0001' u 1:10 t 'alpha2', 'prcc.0001' u 1:12 t 'alpha4', 'prcc.0001' u 1:14 t 'alpha3', 'prcc.0001' u 1:16 t 'alpha10', 'prcc.0001' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: T0"
plot 'prcc.0002' u 1:10 t 'alpha2', 'prcc.0002' u 1:12 t 'alpha4', 'prcc.0002' u 1:14 t 'alpha3', 'prcc.0002' u 1:16 t 'alpha10', 'prcc.0002' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: T1"
plot 'prcc.0003' u 1:10 t 'alpha2', 'prcc.0003' u 1:12 t 'alpha4', 'prcc.0003' u 1:14 t 'alpha3', 'prcc.0003' u 1:16 t 'alpha10', 'prcc.0003' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: T2"
plot 'prcc.0004' u 1:10 t 'alpha2', 'prcc.0004' u 1:12 t 'alpha4', 'prcc.0004' u 1:14 t 'alpha3', 'prcc.0004' u 1:16 t 'alpha10', 'prcc.0004' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: T8"
plot 'prcc.0005' u 1:10 t 'alpha2', 'prcc.0005' u 1:12 t 'alpha4', 'prcc.0005' u 1:14 t 'alpha3', 'prcc.0005' u 1:16 t 'alpha10', 'prcc.0005' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: MR"
plot 'prcc.0006' u 1:10 t 'alpha2', 'prcc.0006' u 1:12 t 'alpha4', 'prcc.0006' u 1:14 t 'alpha3', 'prcc.0006' u 1:16 t 'alpha10', 'prcc.0006' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: MI"
plot 'prcc.0007' u 1:10 t 'alpha2', 'prcc.0007' u 1:12 t 'alpha4', 'prcc.0007' u 1:14 t 'alpha3', 'prcc.0007' u 1:16 t 'alpha10', 'prcc.0007' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: MA"
plot 'prcc.0008' u 1:10 t 'alpha2', 'prcc.0008' u 1:12 t 'alpha4', 'prcc.0008' u 1:14 t 'alpha3', 'prcc.0008' u 1:16 t 'alpha10', 'prcc.0008' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: IG"
plot 'prcc.0009' u 1:10 t 'alpha2', 'prcc.0009' u 1:12 t 'alpha4', 'prcc.0009' u 1:14 t 'alpha3', 'prcc.0009' u 1:16 t 'alpha10', 'prcc.0009' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: I12"
plot 'prcc.0010' u 1:10 t 'alpha2', 'prcc.0010' u 1:12 t 'alpha4', 'prcc.0010' u 1:14 t 'alpha3', 'prcc.0010' u 1:16 t 'alpha10', 'prcc.0010' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: I4"
plot 'prcc.0011' u 1:10 t 'alpha2', 'prcc.0011' u 1:12 t 'alpha4', 'prcc.0011' u 1:14 t 'alpha3', 'prcc.0011' u 1:16 t 'alpha10', 'prcc.0011' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: I10"
plot 'prcc.0012' u 1:10 t 'alpha2', 'prcc.0012' u 1:12 t 'alpha4', 'prcc.0012' u 1:14 t 'alpha3', 'prcc.0012' u 1:16 t 'alpha10', 'prcc.0012' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: FA"
plot 'prcc.0013' u 1:10 t 'alpha2', 'prcc.0013' u 1:12 t 'alpha4', 'prcc.0013' u 1:14 t 'alpha3', 'prcc.0013' u 1:16 t 'alpha10', 'prcc.0013' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: BI"
plot 'prcc.0014' u 1:10 t 'alpha2', 'prcc.0014' u 1:12 t 'alpha4', 'prcc.0014' u 1:14 t 'alpha3', 'prcc.0014' u 1:16 t 'alpha10', 'prcc.0014' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: TB"
plot 'prcc.0015' u 1:10 t 'alpha2', 'prcc.0015' u 1:12 t 'alpha4', 'prcc.0015' u 1:14 t 'alpha3', 'prcc.0015' u 1:16 t 'alpha10', 'prcc.0015' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: T"
plot 'prcc.0016' u 1:10 t 'alpha2', 'prcc.0016' u 1:12 t 'alpha4', 'prcc.0016' u 1:14 t 'alpha3', 'prcc.0016' u 1:16 t 'alpha10', 'prcc.0016' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: M"
plot 'prcc.0017' u 1:10 t 'alpha2', 'prcc.0017' u 1:12 t 'alpha4', 'prcc.0017' u 1:14 t 'alpha3', 'prcc.0017' u 1:16 t 'alpha10', 'prcc.0017' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: MIp"
plot 'prcc.0018' u 1:10 t 'alpha2', 'prcc.0018' u 1:12 t 'alpha4', 'prcc.0018' u 1:14 t 'alpha3', 'prcc.0018' u 1:16 t 'alpha10', 'prcc.0018' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: MITIratio"
plot 'prcc.0019' u 1:10 t 'alpha2', 'prcc.0019' u 1:12 t 'alpha4', 'prcc.0019' u 1:14 t 'alpha3', 'prcc.0019' u 1:16 t 'alpha10', 'prcc.0019' u 1:18 t 'k18', 0 t '' w l, 0.01 t '' w l 
exit
