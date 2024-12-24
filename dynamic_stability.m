clc;
clear;
close all;

a = -.2;
e = -.1;
mo = 20;
R_2 = .24;
sigma = .4;
xteta = e-a;
%%first scenario
zeta1 = 0;
zetateta1 = 0;
%equation
x = 0;
for v1 = .1 : .1 : 4
 x = x+1;   
V1(1,x) = v1;
B1 = [1 xteta (2*zeta1*sigma)/v1 0; xteta R_2 0 (2*zetateta1*R_2)/v1; 0 0 1 0; 0 0 0 1] ;
A1 = -[0 0 sigma^2/v1^2 2/mo; 0 0 0 (R_2/v1^2)-((2/mo)*(.5+a)); -1 0 0 0 ; 0 -1 0 0 ];
[Q,D] = eig(A1,B1);
r1_1(x) = D(1,1);
r2_1(x) = D(2,2);
r3_1(x) = D(3,3); 
r4_1(x) = D(4,4);

end
figure(1);
plot(V1,V1.*imag(r1_1),'ro', 'markerSize',8);
hold on
plot(V1,V1.*imag(r2_1),'bd', 'markerSize',8);
plot(V1,V1.*imag(r3_1),'gs', 'markerSize',8);
plot(V1,V1.*imag(r4_1),'k*', 'markerSize',8);
xlabel('V');
ylabel('\Omega/\Omega_\theta')

grid on
hold off

figure(2);
plot(V1,V1.*real(r1_1),'ro', 'markerSize',8);
hold on
plot(V1,V1.*real(r2_1),'bd', 'markerSize',8);
plot(V1,V1.*real(r3_1),'gs', 'markerSize',8);
plot(V1,V1.*real(r4_1),'k*', 'markerSize',8);
xlabel('V')
ylabel('\Gamma/\Omega_\theta')
grid on
hold off


%%second scenario
zeta2 = 0.05;
zetateta2 = 0;
%equation
y = 0;
for v = .1 : .1 : 4
 y = y+1;   
V(1,y) = v;
B2 = [1 xteta (2*zeta2*sigma)/v 0; xteta R_2 0 (2*zetateta2*R_2)/v; 0 0 1 0; 0 0 0 1]; 
A2 = -[0 0 sigma^2/v^2 2/mo; 0 0 0 (R_2/v^2)-((2/mo)*(.5+a)); -1 0 0 0 ; 0 -1 0 0 ];
[Q2,D2] = eig(A2,B2);
r1_2(y) = D2(1,1);
r2_2(y) = D2(2,2);
r3_2(y) = D2(3,3); 
r4_2(y) = D2(4,4);

end
figure(3);
plot(V,V.*imag(r1_2),'ro', 'markerSize',8);
hold on
plot(V,V.*imag(r2_2),'bd', 'markerSize',8);
plot(V,V.*imag(r3_2),'gs', 'markerSize',8);
plot(V,V.*imag(r4_2),'k*', 'markerSize',8);
xlabel('V');
ylabel('\Omega/\Omega_\theta')

grid on
hold off

figure(4);
plot(V,V.*real(r1_2),'ro', 'markerSize',8);
hold on
plot(V,V.*real(r2_2),'bd', 'markerSize',8);
plot(V,V.*real(r3_2),'gs', 'markerSize',8);
plot(V,V.*real(r4_2),'k*', 'markerSize',8);
xlabel('V')
ylabel('\Gamma/\Omega_\theta')
grid on
hold off


% % third scenario
zeta3 = 0;
zetateta3 = 0.05;
%equation
z = 0;
for v = .1 : .1 : 4
 z = z+1;   
V(1,z) = v;
B3 = [1 xteta (2*zeta3*sigma)/v 0; xteta R_2 0 (2*zetateta3*R_2)/v; 0 0 1 0; 0 0 0 1] ;
A3 = -[0 0 sigma^2/v^2 2/mo; 0 0 0 (R_2/v^2)-((2/mo)*(.5+a)); -1 0 0 0 ; 0 -1 0 0 ];
[Q3,D3] = eig(A3,B3);
r1_3(z) = D3(1,1);
r2_3(z) = D3(2,2);
r3_3(z) = D3(3,3); 
r4_3(z) = D3(4,4);

end
figure(5);
plot(V,V.*imag(r1_3),'ro', 'markerSize',8);
hold on
plot(V,V.*imag(r2_3),'bd', 'markerSize',8);
plot(V,V.*imag(r3_3),'gs', 'markerSize',8);
plot(V,V.*imag(r4_3),'k*', 'markerSize',8);
xlabel('V');
ylabel('\Omega/\Omega_\theta')

grid on
hold off

figure(6);
plot(V,V.*real(r1_3),'ro', 'markerSize',8);
hold on
plot(V,V.*real(r2_3),'bd', 'markerSize',8);
plot(V,V.*real(r3_3),'gs', 'markerSize',8);
plot(V,V.*real(r4_3),'k*', 'markerSize',8);
xlabel('V')
ylabel('\Gamma/\Omega_\theta')
grid on
hold off


% % forth scenario
zeta4 = .05;
zetateta4 = .05;
%equation
o = 0;
for v = .1 : .1 : 4
 o = o+1;   
V(1,o) = v;
B4 = [1 xteta (2*zeta4*sigma)/v 0; xteta R_2 0 (2*zetateta4*R_2)/v; 0 0 1 0; 0 0 0 1] ;
A4 = -[0 0 sigma^2/v^2 2/mo; 0 0 0 (R_2/v^2)-((2/mo)*(.5+a)); -1 0 0 0 ; 0 -1 0 0 ];
[Q4,D4] = eig(A4,B4);
r1_4(o) = D4(1,1);
r2_4(o) = D4(2,2);
r3_4(o) = D4(3,3); 
r4_4(o) = D4(4,4);

end
figure(7);
plot(V,V.*imag(r1_4),'ro', 'markerSize',8);
hold on
plot(V,V.*imag(r2_4),'bd', 'markerSize',8);
plot(V,V.*imag(r3_4),'gs', 'markerSize',8);
plot(V,V.*imag(r4_4),'k*', 'markerSize',8);
xlabel('V');
ylabel('\Omega/\Omega_\theta')

grid on
hold off

figure(8);
plot(V,V.*real(r1_4),'ro', 'markerSize',8);
hold on
plot(V,V.*real(r2_4),'bd', 'markerSize',8);
plot(V,V.*real(r3_4),'gs', 'markerSize',8);
plot(V,V.*real(r4_4),'k*', 'markerSize',8);
xlabel('V')
ylabel('\Gamma/\Omega_\theta')
grid on
hold off
