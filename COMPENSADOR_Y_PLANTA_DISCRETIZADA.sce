clc;
clear;
//Parametros
M=0.696; //masa del carro
m=0.017; //masa del pendulo
l=0.3;    //longitud del pendulo
I=0.0011; //inercia del pendulo
g=9.8;    //aceleracion de gravedad
b=0.001;  //coeficiente de friccion

ap=[0 1 0 0; 
   0 -b*(I+m*l^2)/((M+m)*I+M*m*l^2) m^2*l^2*g/((M+m)*I+M*m*l^2) 0;
   0 0 0 1;
   0 -m*l*b/((M+m)*I+M*m*l^2) g*m*l*(M+m)/((M+m)*I+M*m*l^2) 0];
bp=[0;
   (I+m*l^2)/((M+m)*I+M*m*l^2);
   0;
   m*l/((M+m)*I+M*m*l^2)];
cp=[1 0 0 0;
   0 0 1 0];
dp=[0;0];

//Controlabilidad
cm = [bp ap*bp (ap ^2)*bp (ap ^3)*bp] ;
rcm = rank (cm)
disp ('rango de controlabilidad')
disp(rcm)
//Observabilidad
om = [cp;cp*ap; cp*(ap ^2);cp*( ap^3) ];
rom = rank (om)
disp ('rango de observabilidad')
disp(rom)

//Ganancia LQR
Q=cp'*cp;
R=0.01;
SYS_1=syslin('c',ap,bp,cp,dp)
[G, X] = lqr(SYS_1, Q, R);
disp('ganancia LQR')
disp(G)
//Ganancia del Observador
Q=cp'*cp;
R=[0.01 0; 0 0.01];
SYS_2=syslin('c',ap,bp,cp,dp)
[H, X] = lqe(SYS_2,Q,R);
disp('ganancia Observador')
disp(H)

//Compensador K(s)
Ak=[ap-bp*G,-bp*G;H*cp,ap-H*cp-bp*G];
Bk=[zeros(4,2);H];
Ck=[cp,zeros(2 ,4)];
Dk=[zeros(2,2)];
SYS_K=syslin('c',Ak,Bk,Ck,Dk);

s=poly(0,'s');
s1=syslin('c',1/s);  //integrador 1/s
FT_K=ss2tf(SYS_K)*s1;
SYS_K=tf2ss(FT_K);

T=0.01 //tiempo de muestreo

SYS_Kz=cls2dls(SYS_K,T);   //compensador discretizado en espacio de Estados
FT_Kz=ss2tf(SYS_Kz);        //convertir en forma de transformada z
