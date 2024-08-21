set term postscript color  10 
set output "prccvals.ps"
set ytics 0.2 
set yrange [-1:1] 
set xrange [0:200] 
set xlabel 'Time(Units)' 
set title "PRCC : Poblacion"
plot 'prcc.0001' u 1:9 t 'tasadenacimientos', 'prcc.0001' u 1:11 t 'tasademuertejovenes', 'prcc.0001' u 1:13 t 'tasademuerteadultos', 'prcc.0001' u 1:15 t 'tasadenuevasminas', 'prcc.0001' u 1:17 t 'tasacierreminas', 'prcc.0001' u 1:19 t 'inicioclosemina', 'prcc.0001' u 1:21 t 'numberofclosemine', 'prcc.0001' u 1:23 t 'numerodenuevasminas', 'prcc.0001' u 1:25 t 'cobreproducido', 'prcc.0001' u 1:27 t 'preciocobre', 'prcc.0001' u 1:29 t 'presupuestoporpersona', 0 t '' w l 
set title "PRCC : Minas"
plot 'prcc.0002' u 1:9 t 'tasadenacimientos', 'prcc.0002' u 1:11 t 'tasademuertejovenes', 'prcc.0002' u 1:13 t 'tasademuerteadultos', 'prcc.0002' u 1:15 t 'tasadenuevasminas', 'prcc.0002' u 1:17 t 'tasacierreminas', 'prcc.0002' u 1:19 t 'inicioclosemina', 'prcc.0002' u 1:21 t 'numberofclosemine', 'prcc.0002' u 1:23 t 'numerodenuevasminas', 'prcc.0002' u 1:25 t 'cobreproducido', 'prcc.0002' u 1:27 t 'preciocobre', 'prcc.0002' u 1:29 t 'presupuestoporpersona', 0 t '' w l 
set title "PRCC : Ingresos"
plot 'prcc.0003' u 1:9 t 'tasadenacimientos', 'prcc.0003' u 1:11 t 'tasademuertejovenes', 'prcc.0003' u 1:13 t 'tasademuerteadultos', 'prcc.0003' u 1:15 t 'tasadenuevasminas', 'prcc.0003' u 1:17 t 'tasacierreminas', 'prcc.0003' u 1:19 t 'inicioclosemina', 'prcc.0003' u 1:21 t 'numberofclosemine', 'prcc.0003' u 1:23 t 'numerodenuevasminas', 'prcc.0003' u 1:25 t 'cobreproducido', 'prcc.0003' u 1:27 t 'preciocobre', 'prcc.0003' u 1:29 t 'presupuestoporpersona', 0 t '' w l 
set title "PRCC : Deudas"
plot 'prcc.0004' u 1:9 t 'tasadenacimientos', 'prcc.0004' u 1:11 t 'tasademuertejovenes', 'prcc.0004' u 1:13 t 'tasademuerteadultos', 'prcc.0004' u 1:15 t 'tasadenuevasminas', 'prcc.0004' u 1:17 t 'tasacierreminas', 'prcc.0004' u 1:19 t 'inicioclosemina', 'prcc.0004' u 1:21 t 'numberofclosemine', 'prcc.0004' u 1:23 t 'numerodenuevasminas', 'prcc.0004' u 1:25 t 'cobreproducido', 'prcc.0004' u 1:27 t 'preciocobre', 'prcc.0004' u 1:29 t 'presupuestoporpersona', 0 t '' w l 
set title "PRCC : Presupuesto"
plot 'prcc.0005' u 1:9 t 'tasadenacimientos', 'prcc.0005' u 1:11 t 'tasademuertejovenes', 'prcc.0005' u 1:13 t 'tasademuerteadultos', 'prcc.0005' u 1:15 t 'tasadenuevasminas', 'prcc.0005' u 1:17 t 'tasacierreminas', 'prcc.0005' u 1:19 t 'inicioclosemina', 'prcc.0005' u 1:21 t 'numberofclosemine', 'prcc.0005' u 1:23 t 'numerodenuevasminas', 'prcc.0005' u 1:25 t 'cobreproducido', 'prcc.0005' u 1:27 t 'preciocobre', 'prcc.0005' u 1:29 t 'presupuestoporpersona', 0 t '' w l 
set title "PRCC : PBIpc"
plot 'prcc.0006' u 1:9 t 'tasadenacimientos', 'prcc.0006' u 1:11 t 'tasademuertejovenes', 'prcc.0006' u 1:13 t 'tasademuerteadultos', 'prcc.0006' u 1:15 t 'tasadenuevasminas', 'prcc.0006' u 1:17 t 'tasacierreminas', 'prcc.0006' u 1:19 t 'inicioclosemina', 'prcc.0006' u 1:21 t 'numberofclosemine', 'prcc.0006' u 1:23 t 'numerodenuevasminas', 'prcc.0006' u 1:25 t 'cobreproducido', 'prcc.0006' u 1:27 t 'preciocobre', 'prcc.0006' u 1:29 t 'presupuestoporpersona', 0 t '' w l 
set title "PRCC : DEUpc"
plot 'prcc.0007' u 1:9 t 'tasadenacimientos', 'prcc.0007' u 1:11 t 'tasademuertejovenes', 'prcc.0007' u 1:13 t 'tasademuerteadultos', 'prcc.0007' u 1:15 t 'tasadenuevasminas', 'prcc.0007' u 1:17 t 'tasacierreminas', 'prcc.0007' u 1:19 t 'inicioclosemina', 'prcc.0007' u 1:21 t 'numberofclosemine', 'prcc.0007' u 1:23 t 'numerodenuevasminas', 'prcc.0007' u 1:25 t 'cobreproducido', 'prcc.0007' u 1:27 t 'preciocobre', 'prcc.0007' u 1:29 t 'presupuestoporpersona', 0 t '' w l 
exit
