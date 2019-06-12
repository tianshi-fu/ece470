clear
clc
%path

S = [0.00 0.00; 0.00 0.00; 0.00 -1.00; 0.00 -2.00; 1.00 0.00; 0.00 0.00];
M = [0.00 0.00 -1.00 0.00; 0.00 -1.00 0.00 0.00; -1.00 0.00 0.00 -4.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 0.00 0.00 0.00; 0.00 0.00 2.00 0.00; 0.00 -2.00 -4.00 -4.00];
r_robot = [0.90 0.90 0.90 0.90];
p_obstacle = [-3.03 4.06 4.22 -1.23 2.39 -4.37 -3.50 3.33 -2.64 -2.58 -2.37 2.50 4.69 -3.55 3.91 -3.22 0.24 -2.98 -3.69 1.38 2.85; -1.78 -1.33 2.90 -0.96 3.88 -1.48 2.38 -2.50 3.90 -3.32 -2.17 -4.72 -1.80 4.51 -3.06 0.24 4.57 -4.92 -0.33 -2.70 2.85; 4.79 2.62 4.20 3.37 4.74 4.55 -1.29 -2.45 0.17 -3.15 -2.35 3.20 -0.03 -3.64 3.01 3.67 1.16 -4.66 -1.10 3.92 -3.49];
r_obstacle = [1.57 2.37 3.93 1.35 3.88 2.74 1.79 1.43 1.25 1.33 0.81 2.42 2.24 2.37 2.63 3.44 1.81 3.82 0.90 3.79 1.29];
theta_a = [0.93; -2.59];
theta_b = [-2.20; -1.85];

S1 = S(:,1);
S2 = S(:,2);


D = norm(theta_b-theta_a);
esp = 0.1;
n = ceil(D/esp)+1;
s = linspace(0,1,n);


for e = 1:n
theta = (1-s(e))*theta_a + s(e) * theta_b;
p1 = p_robot(:,1);
p2 = p_robot(:,2);
p3 = expm(ssV(S1,theta(1)))*[p_robot(:,3);1];
p4 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*[p_robot(:,4);1];

p = [p1 p2 p3(1:3) p4(1:3)];

    
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
        c(e) = s(e);
    else
        c(e) = 0;
end

end
disp(mat2str(c,6))


