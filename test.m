clear
clc
rng('shuffle')
theta = randn(7,1);

T_1in0 = [-0.91699320 -0.37170461 0.14477280 -0.44194233; 0.05928019 0.23191640 0.97092772 -3.16216328; -0.39447350 0.89891627 -0.19063104 6.12809008; 0.00000000 0.00000000 0.00000000 1.00000000];
M = [0.00000000 -1.00000000 0.00000000 0.00000000; 0.00000000 0.00000000 1.00000000 -2.00000000; -1.00000000 0.00000000 0.00000000 6.00000000; 0.00000000 0.00000000 0.00000000 1.00000000];
S = [0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000; 0.00000000 0.00000000 0.00000000 -1.00000000 -1.00000000 0.00000000 0.00000000; 0.00000000 0.00000000 1.00000000 0.00000000 0.00000000 0.00000000 -1.00000000; 0.00000000 0.00000000 4.00000000 6.00000000 6.00000000 1.00000000 2.00000000; 1.00000000 1.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000; 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000];
S1 = S(:,1);
S2 = S(:,2);
S3 = S(:,3);
S4 = S(:,4);
S5 = S(:,5);
S6 = S(:,6);
S7 = S(:,7);

n = 2;

while n > 0.01
    TC = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5)))*expm(ssV(S6,theta(6)))*expm(ssV(S7,theta(7)))*M;
    ssVs = logm(T_1in0*(TC^-1));
    Vs = [ssVs(7); -ssVs(3); ssVs(2); ssVs(1:3,4)];
    n = norm(Vs);
    J1 = S1;
    J2 = adjo(expm(ssV(S1,theta(1))))*S2;
    J3 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2))))*S3;
    J4 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3))))*S4;
    J5 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4))))*S5;
    J6 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5))))*S6;
    J7 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5)))*expm(ssV(S6,theta(6))))*S7;
    J = [J1 J2 J3 J4 J5 J6 J7];
    thetadot = inv(J' * J + 0.1 * eye(7))*J'*Vs;
    theta = theta + thetadot;
end

disp(mat2str(theta,4))