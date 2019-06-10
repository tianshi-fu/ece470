function A = IK_final(next_T,M,old_theta)
rng('shuffle')
% M =[1 0 0 -1.3056; 0 1 0 0.4; 0 0 1 0.96915; 0 0 0 1];
q1 = [-1;0.4;0.4225];
q2 = [-1.1117;0.39997;0.4269];
q3 = [-1.1117;0.40004;0.67055];
q4 = [-1.1117;0.4;0.8838];
q5 = [-1.1123;0.4;0.96802];
q6 = [-1.1117;0.4;0.96915];
a1 = [0;0;1];
a2 = [-1;0;0];
a3 = [-1;0;0];
a4 = [-1;0;0];
a5 = [0;0;1];
a6 = [-1;0;0];
S1 = [a1;-ss31(a1)*q1];
S2 = [a2;-ss31(a2)*q2];
S3 = [a3;-ss31(a3)*q3];
S4 = [a4;-ss31(a4)*q4];
S5 = [a5;-ss31(a5)*q5];
S6 = [a6;-ss31(a6)*q6];
theta = old_theta;

n = 2;
tic; %start timer
e = 0; %time elapsed 
 
while n > 0.01 && e < 5 %run algorithm only if norm of body twist is larger than 0.01 and elapsed time is lesser than 4 seconds
    TC = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5)))*expm(ssV(S6,theta(6)))*M;
    ssVs = logm(next_T*(TC^-1));
    Vs = [ssVs(7); -ssVs(3); ssVs(2); ssVs(1:3,4)];
    n = norm(Vs);
    J1 = S1;
    J2 = adjo(expm(ssV(S1,theta(1))))*S2;
    J3 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2))))*S3;
    J4 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3))))*S4;
    J5 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4))))*S5;
    J6 = adjo(expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5))))*S6;
    J = [J1 J2 J3 J4 J5 J6];
    thetadot = inv(J' * J + 0.01 * eye(6))*J'*Vs;
    theta = theta + thetadot;
    
    e = toc;
    
    mesg = ['Processing (',num2str(e),' seconds)'];
    disp(mesg)
   
end


if e < 4 % pass back joint variables only if algorithm is completed under 4 seconds
  A = theta;
else
  disp('>>>>Impossible Pose!!<<<<')
  A = []; % pass back empty arrary if elapsed time is more than 4 seconds
end

end