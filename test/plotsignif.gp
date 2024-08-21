set term postscript color  10 
set output "signif.ps"
set yrange [0:0.3] 
set ytics 0.01 
set xrange [0:200] 
set xlabel 'Time(Units)' 
set title "PRCC Significance: Poblacion"
plot 'prcc.0001' u 1:10 t 'tasadenacimientos', 'prcc.0001' u 1:12 t 'tasademuertejovenes', 'prcc.0001' u 1:14 t 'tasademuerteadultos', 'prcc.0001' u 1:16 t 'tasadenuevasminas', 'prcc.0001' u 1:18 t 'tasacierreminas', 'prcc.0001' u 1:20 t 'inicioclosemina', 'prcc.0001' u 1:22 t 'numberofclosemine', 'prcc.0001' u 1:24 t 'numerodenuevasminas', 'prcc.0001' u 1:26 t 'cobreproducido', 'prcc.0001' u 1:28 t 'preciocobre', 'prcc.0001' u 1:30 t 'presupuestoporpersona', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: Minas"
plot 'prcc.0002' u 1:10 t 'tasadenacimientos', 'prcc.0002' u 1:12 t 'tasademuertejovenes', 'prcc.0002' u 1:14 t 'tasademuerteadultos', 'prcc.0002' u 1:16 t 'tasadenuevasminas', 'prcc.0002' u 1:18 t 'tasacierreminas', 'prcc.0002' u 1:20 t 'inicioclosemina', 'prcc.0002' u 1:22 t 'numberofclosemine', 'prcc.0002' u 1:24 t 'numerodenuevasminas', 'prcc.0002' u 1:26 t 'cobreproducido', 'prcc.0002' u 1:28 t 'preciocobre', 'prcc.0002' u 1:30 t 'presupuestoporpersona', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: Ingresos"
plot 'prcc.0003' u 1:10 t 'tasadenacimientos', 'prcc.0003' u 1:12 t 'tasademuertejovenes', 'prcc.0003' u 1:14 t 'tasademuerteadultos', 'prcc.0003' u 1:16 t 'tasadenuevasminas', 'prcc.0003' u 1:18 t 'tasacierreminas', 'prcc.0003' u 1:20 t 'inicioclosemina', 'prcc.0003' u 1:22 t 'numberofclosemine', 'prcc.0003' u 1:24 t 'numerodenuevasminas', 'prcc.0003' u 1:26 t 'cobreproducido', 'prcc.0003' u 1:28 t 'preciocobre', 'prcc.0003' u 1:30 t 'presupuestoporpersona', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: Deudas"
plot 'prcc.0004' u 1:10 t 'tasadenacimientos', 'prcc.0004' u 1:12 t 'tasademuertejovenes', 'prcc.0004' u 1:14 t 'tasademuerteadultos', 'prcc.0004' u 1:16 t 'tasadenuevasminas', 'prcc.0004' u 1:18 t 'tasacierreminas', 'prcc.0004' u 1:20 t 'inicioclosemina', 'prcc.0004' u 1:22 t 'numberofclosemine', 'prcc.0004' u 1:24 t 'numerodenuevasminas', 'prcc.0004' u 1:26 t 'cobreproducido', 'prcc.0004' u 1:28 t 'preciocobre', 'prcc.0004' u 1:30 t 'presupuestoporpersona', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: Presupuesto"
plot 'prcc.0005' u 1:10 t 'tasadenacimientos', 'prcc.0005' u 1:12 t 'tasademuertejovenes', 'prcc.0005' u 1:14 t 'tasademuerteadultos', 'prcc.0005' u 1:16 t 'tasadenuevasminas', 'prcc.0005' u 1:18 t 'tasacierreminas', 'prcc.0005' u 1:20 t 'inicioclosemina', 'prcc.0005' u 1:22 t 'numberofclosemine', 'prcc.0005' u 1:24 t 'numerodenuevasminas', 'prcc.0005' u 1:26 t 'cobreproducido', 'prcc.0005' u 1:28 t 'preciocobre', 'prcc.0005' u 1:30 t 'presupuestoporpersona', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: PBIpc"
plot 'prcc.0006' u 1:10 t 'tasadenacimientos', 'prcc.0006' u 1:12 t 'tasademuertejovenes', 'prcc.0006' u 1:14 t 'tasademuerteadultos', 'prcc.0006' u 1:16 t 'tasadenuevasminas', 'prcc.0006' u 1:18 t 'tasacierreminas', 'prcc.0006' u 1:20 t 'inicioclosemina', 'prcc.0006' u 1:22 t 'numberofclosemine', 'prcc.0006' u 1:24 t 'numerodenuevasminas', 'prcc.0006' u 1:26 t 'cobreproducido', 'prcc.0006' u 1:28 t 'preciocobre', 'prcc.0006' u 1:30 t 'presupuestoporpersona', 0 t '' w l, 0.01 t '' w l 
set title "PRCC Significance: DEUpc"
plot 'prcc.0007' u 1:10 t 'tasadenacimientos', 'prcc.0007' u 1:12 t 'tasademuertejovenes', 'prcc.0007' u 1:14 t 'tasademuerteadultos', 'prcc.0007' u 1:16 t 'tasadenuevasminas', 'prcc.0007' u 1:18 t 'tasacierreminas', 'prcc.0007' u 1:20 t 'inicioclosemina', 'prcc.0007' u 1:22 t 'numberofclosemine', 'prcc.0007' u 1:24 t 'numerodenuevasminas', 'prcc.0007' u 1:26 t 'cobreproducido', 'prcc.0007' u 1:28 t 'preciocobre', 'prcc.0007' u 1:30 t 'presupuestoporpersona', 0 t '' w l, 0.01 t '' w l 
exit
