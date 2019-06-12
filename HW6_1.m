clear
clc
%multi-path

S = [0.00 0.00 0.00 0.00; 0.00 0.00 0.00 0.00; 0.00 0.00 1.00 0.00; 1.00 0.00 0.00 -1.00; 0.00 -1.00 -8.00 0.00; 0.00 0.00 0.00 0.00];
M = [-1.00 0.00 0.00 4.00; 0.00 1.00 0.00 0.00; 0.00 0.00 -1.00 0.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 2.00 8.00 8.00 6.00 4.00; 0.00 0.00 0.00 0.00 0.00 0.00; 0.00 -2.00 -2.00 0.00 0.00 0.00];
r_robot = [0.90 0.90 0.90 0.90 0.90 0.90];
p_obstacle = [-2.31 -3.53 -2.70 -4.97 1.59 3.06 -0.60 -1.67 2.75 2.54 0.90 1.41 -4.24 3.09; 0.23 4.08 4.63 2.36 3.51 0.12 -1.65 2.84 -3.30 2.17 -4.78 3.04 2.83 3.26; 4.69 -0.67 4.62 1.36 -2.78 4.82 4.45 0.74 -4.59 -4.43 1.19 -4.40 -2.32 -0.61];
r_obstacle = [3.64 3.83 0.65 4.57 1.87 1.41 2.80 2.44 3.03 1.34 0.71 1.71 1.20 1.18];
theta_start = [2.80 -1.07 2.84 -3.04 1.30 2.90 -2.10 0.76 1.05; -2.37 -0.38 -0.60 2.89 -1.56 -0.53 1.62 1.48 0.40; -2.66 0.12 2.26 1.60 0.43 -1.38 1.01 -0.65 2.69; 1.73 0.25 0.87 1.48 2.96 2.45 0.13 1.15 0.67];
theta_goal = [-1.47 0.73 2.96 2.65 2.49 2.67 -2.17 1.68 1.03; -2.65 -1.82 -0.16 1.63 3.14 1.86 1.23 -2.09 -1.67; 3.01 -0.70 0.70 -1.42 -1.49 -2.28 1.30 2.30 -0.86; 1.49 0.66 2.44 0.73 1.83 0.17 1.86 3.06 0.78];

[~,npaths]=size(theta_start);
S1 = S(:,1);
S2 = S(:,2);
S3 = S(:,3);
S4 = S(:,4);
% S5 = S(:,5);
% S6 = S(:,6);

for x = 1:npaths
D = norm(theta_goal(:,x)-theta_start(:,x));
esp = 0.1; %resolution
n = ceil(D/esp)+1;
s = linspace(0,1,n);


for e = 1:n
theta = (1-s(e))*theta_start(:,x) + s(e) * theta_goal(:,x);
p1 = p_robot(:,1);
p2 = p_robot(:,2);
p3 = expm(ssV(S1,theta(1)))*[p_robot(:,3);1];
p4 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*[p_robot(:,4);1];
p5 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*[p_robot(:,5);1];
p6 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*[p_robot(:,6);1];
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
        c(x,e) = s(e);
    else
        c(x,e) = 0;
end

end
end
for i=1:npaths 
A(i) = max(c(i,:));
end
disp(mat2str(A,6))