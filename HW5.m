clear
clc
%nconfigs

S = [0.00 0.00 -1.00 0.00 0.00; -1.00 0.00 0.00 1.00 1.00; 0.00 0.00 0.00 0.00 0.00; -2.00 0.00 0.00 4.00 6.00; 0.00 0.00 4.00 0.00 0.00; 0.00 -1.00 -2.00 0.00 0.00];
M = [0.00 0.00 1.00 0.00; -1.00 0.00 0.00 0.00; 0.00 -1.00 0.00 -8.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 0.00 0.00 0.00 0.00 0.00 0.00; 0.00 0.00 -2.00 -2.00 0.00 0.00 0.00; 0.00 -2.00 -2.00 -4.00 -4.00 -6.00 -8.00];
r_robot = [0.90 0.90 0.90 0.90 0.90 0.90 0.90];
p_obstacle = [-4.36 -7.74 2.24 -6.65 0.62 3.05 8.32 4.98 -1.20 9.72 5.10 -7.17 -8.59 7.98 -0.56 -1.49; 0.15 5.88 -9.81 3.25 0.62 3.59 7.25 -9.38 3.87 5.09 -2.47 -2.14 -0.68 -6.47 6.47 8.31; 2.32 -3.03 -2.50 -0.74 2.74 9.66 1.09 -8.78 -1.34 -5.05 4.27 8.92 -5.64 -2.29 -6.39 1.59];
r_obstacle = [1.97 3.51 1.43 3.81 0.83 3.41 3.36 2.46 1.67 2.59 3.27 4.56 4.17 1.64 0.62 4.70];
theta = [2.02 1.63 -2.67 -0.32 2.62 1.04 1.22 2.36 -2.52 2.29 2.18; 1.07 2.05 1.55 1.54 1.96 0.77 -1.36 -0.02 2.73 0.66 -2.08; 0.84 2.61 0.32 2.07 1.40 1.46 -1.14 -0.12 -3.13 -0.85 0.23; 1.81 1.62 0.06 1.80 2.60 2.18 -0.18 -1.84 -0.55 -2.52 -2.92; -2.02 1.66 0.80 0.60 -1.85 -1.19 2.71 -1.20 2.19 -2.07 -2.89];

[~,nconfigs]=size(theta);
S1 = S(:,1);
S2 = S(:,2);
S3 = S(:,3);
S4 = S(:,4);
% S5 = S(:,5);
% S6 = S(:,6);

for e = 1:nconfigs
p1 = p_robot(:,1);
p2 = p_robot(:,2);
p3 = expm(ssV(S1,theta(1,e)))*[p_robot(:,3);1];
p4 = expm(ssV(S1,theta(1,e)))*expm(ssV(S2,theta(2,e)))*[p_robot(:,4);1];
p5 = expm(ssV(S1,theta(1,e)))*expm(ssV(S2,theta(2,e)))*expm(ssV(S3,theta(3,e)))*[p_robot(:,5);1];
p6 = expm(ssV(S1,theta(1,e)))*expm(ssV(S2,theta(2,e)))*expm(ssV(S3,theta(3,e)))*expm(ssV(S4,theta(4,e)))*[p_robot(:,6);1];
% p7 = expm(ssV(S1,theta(1,e)))*expm(ssV(S2,theta(2,e)))*expm(ssV(S3,theta(3,e)))*expm(ssV(S4,theta(4,e)))*expm(ssV(S5,theta(5,e)))*[p_robot(:,7);1];
% p8 = expm(ssV(S1,theta(1,e)))*expm(ssV(S2,theta(2,e)))*expm(ssV(S3,theta(3,e)))*expm(ssV(S4,theta(4,e)))*expm(ssV(S5,theta(5,e)))*expm(ssV(S6,theta(6,e)))*[p_robot(:,8);1];
p = [p1 p2 p3(1:3) p4(1:3) p5(1:3) p6(1:3)];

    
[~,nspheres]=size(p_robot);
[~,nobst]=size(p_obstacle);
combo = combnk(1:nspheres,2);
[ncombo,~]=size(combo);
r = mean(r_robot);

for f = 1:ncombo %for self-collision
d_s(f) = norm(p(:,combo(f,1))-p(:,combo(f,2)))-r-r;    
end

for v = 1:nspheres %for obstacle collision
for g = 1:nobst 
d_o(v,g) = norm(p(:,v)-p_obstacle(:,g))-r-r_obstacle(g);
end
end

if min(d_s) <= 0 || min(d_o(:)) <= 0
        c(e) = 1;
    else
        c(e) = 0;
end

end
disp(mat2str(c))


