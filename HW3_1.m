clear
clc

S = [0.00 0.00 0.00 0.00 0.00; -1.00 0.00 0.00 0.00 0.00; 0.00 0.00 0.00 0.00 0.00; 0.00 0.00 0.00 0.00 1.00; 0.00 -1.00 0.00 0.00 0.00; 0.00 0.00 1.00 1.00 0.00];
M = [0.00 0.00 1.00 -4.00; 0.00 1.00 0.00 0.00; -1.00 0.00 0.00 0.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 0.00 -2.00 -2.00 -6.00 -6.00 -4.00; 0.00 -2.00 -2.00 -4.00 -2.00 0.00 0.00; 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
r_robot = [0.90 0.90 0.90 0.90 0.90 0.90 0.90];
p_obstacle = [8.14 5.96 5.01 5.83 8.36 9.38 7.26 5.48 8.86 1.92 -2.91 -6.87 -7.88 6.91 -4.60 1.60 8.73 3.43 -1.87; 1.40 5.75 6.85 -6.06 0.27 -1.78 5.70 -4.11 -9.63 8.01 8.77 4.09 4.76 7.49 -3.26 2.27 -2.56 -0.23 7.70; -3.79 -0.38 -7.91 -2.20 0.29 2.46 2.78 -4.27 9.01 -6.43 1.12 -3.33 -6.02 7.71 3.83 7.44 7.18 -3.31 1.00];
r_obstacle = [1.45 3.24 4.28 1.27 0.56 0.96 3.44 2.17 4.68 3.60 3.81 0.71 2.66 4.99 2.37 1.34 4.19 0.62 3.66];
theta = [-2.39; -1.92; -0.25; -1.98; 1.82];


S1 = S(:,1);
S2 = S(:,2);
S3 = S(:,3);
S4 = S(:,4);
S5 = S(:,5);
% S6 = S(:,6);

[~,y]=size(p_robot);
[~,u]=size(p_obstacle);
[~,x]=size(theta);
c = zeros(1,x);
n = [];
m = [];
r = mean(r_robot);

for e = 1:x
thetac = theta(:,e);
p1 = [0;0;0];
p2 = p_robot(:,2);
p3 = expm(ssV(S1,thetac(1)))*[p_robot(:,3);1];
p4 = expm(ssV(S1,thetac(1)))*expm(ssV(S2,thetac(2)))*[p_robot(:,4);1];
p5 = expm(ssV(S1,thetac(1)))*expm(ssV(S2,thetac(2)))*expm(ssV(S3,thetac(3)))*[p_robot(:,5);1];
p6 = expm(ssV(S1,thetac(1)))*expm(ssV(S2,thetac(2)))*expm(ssV(S3,thetac(3)))*expm(ssV(S4,thetac(4)))*[p_robot(:,6);1];
p7 = expm(ssV(S1,thetac(1)))*expm(ssV(S2,thetac(2)))*expm(ssV(S3,thetac(3)))*expm(ssV(S4,thetac(4)))*expm(ssV(S5,thetac(5)))*[p_robot(:,7);1];
% p8 = expm(ssV(S1,thetac(1)))*expm(ssV(S2,thetac(2)))*expm(ssV(S3,thetac(3)))*expm(ssV(S4,thetac(4)))*expm(ssV(S5,thetac(5)))*expm(ssV(S6,thetac(6)))*M(:,4);

p = [p1 p2 p3(1:3) p4(1:3) p5(1:3) p6(1:3) p7(1:3)];

combo = combnk(1:y,2);
[z,~]=size(combo);
for f = 1:z %norm for self-collision
n(e,f) = norm(p(:,combo(f,1))-p(:,combo(f,2)));    
end   

for v = 1:y
for g = 1:u*y %norm for obstacle collision
m(v,g) = norm(p(:,g)-p_obstacle(:,g));    
end
end

min_n = min(n(e,:));
[min_m,ind] = min(m(:));
d_s = min_n-r-r;
d_o = min_m-r-r_obstacle(ind);
d = [d_s d_o];
if all(d) < 0
        c(e) = 1;
    else
        c(e) = 0;
end

end
disp(mat2str(c,4))