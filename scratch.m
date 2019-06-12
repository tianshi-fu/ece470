clear
clc

q_in1 = [0.58691733 0.92153560 0.65097551 1.61176785 -1.18396295 -0.09655846 -1.33698044 1.22678897 0.55164859; -1.06410912 0.03621220 0.35813227 0.74189181 -0.16490916 -1.08283899 0.91418564 0.73693104 1.35440634; -3.46257519 -0.45107053 1.14296809 -0.29753986 0.03899091 0.58205288 -0.48497808 0.48099701 -0.64936159; 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000];
q_in0 = [4.22936442 1.77535809 0.62787083 1.20648079 1.37238522 2.01818683 0.72088699 0.74420253 0.72012221; -1.88884567 0.12282995 1.29220321 -0.06621621 0.49784252 1.51387305 -0.52084253 0.57181570 -0.75851685; -3.02574603 -3.59394571 -3.41055108 -4.39308030 -1.49377605 -2.42682602 -1.50908163 -4.02761196 -3.44072501; 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000];

rng shuffle
R = rand(3,3);
p = rand(3,1);
M = [R p; 0 0 0 1];



n = 2;
[~,npoints]=size(q_in0);
sum_A = zeros(6);
sum_B = zeros(6,1);
ident = eye(6);

while n > 1e-8
      b = q_in0(:,1)- M*q_in1(:,1); 
             
      J1 = ss61(ident(:,1))*M*q_in1(:,1);
      J2 = ss61(ident(:,2))*M*q_in1(:,1);
      J3 = ss61(ident(:,3))*M*q_in1(:,1);
      J4 = ss61(ident(:,4))*M*q_in1(:,1);
      J5 = ss61(ident(:,5))*M*q_in1(:,1);
      J6 = ss61(ident(:,6))*M*q_in1(:,1);
       
      J = [J1 J2 J3 J4 J5 J6];
   
   
             
   u = inv(J' * J + 0.1 * eye(6)) * J' * b;
   
   
   
   M = expm(ss61(u))*M;           
   n=norm(u);
   
    disp(n)
    
end
disp(mat2str(M,5));