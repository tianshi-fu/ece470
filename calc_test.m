function [A,B,C]=calc_test(M,K,PI_0,q,p)
ident = eye(6);
r = K*PI_0*M*p;
eta = [r(1)/r(3); r(2)/r(3); 1];

C=norm(q - eta)^2;

b = q - eta;

deri_eta = [1/r(3) 0 -r(1)/(r(3)^2); 0 1/r(3) -r(2)/(r(3)^2); 0 0 0];

J1 = ss61(ident(:,1))*M*p;
J2 = ss61(ident(:,2))*M*p;
J3 = ss61(ident(:,3))*M*p;
J4 = ss61(ident(:,4))*M*p;
J5 = ss61(ident(:,5))*M*p;
J6 = ss61(ident(:,6))*M*p;
J = deri_eta*K*PI_0*[J1 J2 J3 J4 J5 J6];


A = J' * J;
B = J' * b;




end