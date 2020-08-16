// Se cargan las matrices ABCD obtenidas del Scilab/////////////////////////////////////

ap=A;
bp=B;
cp=C;
dp=D;
// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap,bp)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap, cp)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a
/* Plot singular values of LTI the model */
G = syslin('c', ap, bp, cp, dp);
tr = trzeros(G)
w = logspace(-3,3);
sv = svplot(G,w);
scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("valores singulares (SCILAB)","Frequency (rad/s)", "Amplitude (dB)");
// Scaling

//graficar polos y ceros de FT SCILAB
FT=ss2tf(sys);
scf(2);
plzr(FT); 
xtitle("Polos y ceros FT SCILAB");
//matrices obtenidas analiticamente:////////////////////////////////////////
exec('edsonjParameters.sce', -1);

ap1=[0 1 0 0; 
   0 -b*(I+m*l^2)/((M+m)*I+M*m*l^2) m^2*l^2*g/((M+m)*I+M*m*l^2) 0;
   0 0 0 1;
   0 -m*l*b/((M+m)*I+M*m*l^2) g*m*l*(M+m)/((M+m)*I+M*m*l^2) 0];
bp1=[0;
   (I+m*l^2)/((M+m)*I+M*m*l^2);
   0;
   m*l/((M+m)*I+M*m*l^2)];
cp1=[1 0 0 0;
   0 0 1 0];
dp1=[0;0];

// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc1 = cont_mat(ap1,bp1)
rankCc1=rank(Cc1)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O1 = obsv_mat(ap1, cp1)
rankO1=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a
/* Plot singular values of LTI the model */
G1 = syslin('c', ap1, bp1, cp1, dp1);
tr1 = trzeros(G1)
w1 = logspace(-3,3);
sv1 = svplot(G1,w1);
scf(3);
plot2d("ln", w1, 20*log(sv1')/log(10))
xgrid(12)
xtitle("valores singulares (ANALITICAMENTE)","Frequency (rad/s)", "Amplitude (dB)");
// Scaling
//graficar polos y ceros de FT ANALITICA
FT1=ss2tf(G1);
scf(4);
plzr(FT1); 
xtitle("Polos y ceros FT ANALITICA");
