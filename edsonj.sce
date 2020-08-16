function xdot=edsonj(u1,u2,u3,u4,u5)
// This is nonlinear EDSON-J model

// Load the parameters
exec('edsonjParameters.sce', -1);

// state variables
x1=u1;		
x2=u2;
x3=u3;
x4=u4;
G=[10 8.6733348  -43.840001  -9.7685989];
F=-G*[x1;x2;x3;x4];

xdot=[x2;
      (-b*(I+m*l^2)*x2+m^2*l^2*g*sin(x3)*cos(x3)-m*l*(I+m*l^2)*x4^2*sin(x3)+(I+m*l^2)*F)/((M+m)*(I+m*l^2)-m^2*l^2*cos(x3)^2);
      x4;
      m*l*cos(x3)*(-b*(I+m*l^2)*x2+m^2*l^2*g*sin(x3)*cos(x3)-m*l*(I+m*l^2)*x4^2*sin(x3)+(I+m*l^2)*F)/((M+m)*(I+m*l^2)^2-(I+m*l^2)*m^2*l^2*cos(x3)^2)+m*g*l*sin(x3)/(I+m*l^2)];
endfunction
