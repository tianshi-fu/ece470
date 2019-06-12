function A  = invK2(T_1in0, M, S)
rng('shuffle')
[~,joints] = size(S);
%modify to suit specific number of joints----------------------------------
B1 = adjo(M^-1)*S(:,1);
B2 = adjo(M^-1)*S(:,2);
B3 = adjo(M^-1)*S(:,3);
B4 = adjo(M^-1)*S(:,4);
B5 = adjo(M^-1)*S(:,5);
B6 = adjo(M^-1)*S(:,6);

theta = randn(joints,1);
TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)));
%--------------------------------------------------------------------------
ssVb = logm((TC^-1)*T_1in0);
Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
n = norm(Vb);

while n > 0.01
    n = norm(Vb);
    %modify to suit specific number of joints----------------------------------
    B1 = adjo(M^-1)*S(:,1);
    B2 = adjo(M^-1)*S(:,2);
    B3 = adjo(M^-1)*S(:,3);
    B4 = adjo(M^-1)*S(:,4);
    B5 = adjo(M^-1)*S(:,5);
    B6 = adjo(M^-1)*S(:,6);
   
    %Jacobian Calculation--------------------------------------------------
    J6 = B6;
    J5 = adjo(expm(-ssV(B6,theta(6))))*B5;
    J4 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
    J3 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
    J2 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
    J1 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;

    J = [J1 J2 J3 J4 J5 J6];
    %----------------------------------------------------------------------
%--------------------------------------------------------------------------
    mu = 0.1;
    thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
    theta = theta + thetadot;
%modify to suit specific number of joints----------------------------------
TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)));
%--------------------------------------------------------------------------
    
    
    ssVb = logm((TC^-1)*T_1in0);
    Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];

end
 A = theta;

end