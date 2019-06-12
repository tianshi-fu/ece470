function collision_state = colli_nj(S,theta_a,theta_b,p_robot,r_robot,p_obstacle,r_obstacle)
S1 = S(:,1);
S2 = S(:,2);
S3 = S(:,3);
S4 = S(:,4);
S5 = S(:,5);
S6 = S(:,6);
S7 = S(:,7);

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
p5 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*[p_robot(:,5);1];
p6 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*[p_robot(:,6);1];
p7 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5)))*[p_robot(:,7);1];
p8 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5)))*expm(ssV(S6,theta(6)))*[p_robot(:,8);1];
p9 = expm(ssV(S1,theta(1)))*expm(ssV(S2,theta(2)))*expm(ssV(S3,theta(3)))*expm(ssV(S4,theta(4)))*expm(ssV(S5,theta(5)))*expm(ssV(S6,theta(6)))*expm(ssV(S7,theta(7)))*[p_robot(:,9);1];
p = [p1 p2 p3(1:3) p4(1:3) p5(1:3) p6(1:3) p7(1:3) p8(1:3) p9(1:3)];

    
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

if 1 && all(c == 0)   
    collision_state=0;
else
    collision_state=1;
end
end
end