function A = invK(T_1in0, M, S)
rng('shuffle')
%modify to suit specific number of joints----------------------------------
B1 = adjo(M^-1)*S(:,1);
B2 = adjo(M^-1)*S(:,2);
B3 = adjo(M^-1)*S(:,3);
B4 = adjo(M^-1)*S(:,4);
B5 = adjo(M^-1)*S(:,5);
B6 = adjo(M^-1)*S(:,6);
B7 = adjo(M^-1)*S(:,7);
B8 = adjo(M^-1)*S(:,8);
B9 = adjo(M^-1)*S(:,9);
theta = randn(9,1);
TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)));
%--------------------------------------------------------------------------
ssVb = logm((TC^-1)*T_1in0);
Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
w = Vb(1:3);
v = Vb(4:6);
j = norm(w);
k = norm(v);

e = 0;
t1 = clock;
 
while (j > 0.01 && k > 0.01 )
while (e < 1)
    w = Vb(1:3);
    v = Vb(4:6);
    j = norm(w);
    k = norm(v);
%modify to suit specific number of joints----------------------------------
    B1 = adjo(M^-1)*S(:,1);
    B2 = adjo(M^-1)*S(:,2);
    B3 = adjo(M^-1)*S(:,3);
    B4 = adjo(M^-1)*S(:,4);
    B5 = adjo(M^-1)*S(:,5);
    B6 = adjo(M^-1)*S(:,6);
    B7 = adjo(M^-1)*S(:,7);
    B8 = adjo(M^-1)*S(:,8);
    B9 = adjo(M^-1)*S(:,9);
    %Jacobian Calculation--------------------------------------------------
    J9 = B9;
    J8 = adjo(expm(-ssV(B9,theta(9))))*B8;
    J7 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8))))*B7;
    J6 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7))))*B6;
    J5 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6))))*B5;
    J4 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
    J3 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
    J2 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
    J1 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
    J = [J1 J2 J3 J4 J5 J6 J7 J8 J9];
    %----------------------------------------------------------------------
%--------------------------------------------------------------------------
    thetadot = pinv(J)*Vb;
    theta = theta + thetadot;
    
    
%modify to suit specific number of joints----------------------------------
    TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)));
%--------------------------------------------------------------------------
    
    
    ssVb = logm((TC^-1)*T_1in0);
    Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
    t2 = clock;
    e = t2 - t1;
end
end

 disp(norm(w))
 disp(norm(v))
 disp(any(theta > 2*pi))
 A = theta;   
end