function A = coli_test(S,p_robot,r_robot,p_obstacle,r_obstacle,theta_start,theta_goal)

S1 = S(:,1);
S2 = S(:,2);
S3 = S(:,3);
S4 = S(:,4);



D = norm(theta_goal-theta_start);
nsamples = ceil(D/0.1)+1;
s = linspace(0,1,nsamples);

for e = 1:nsamples
theta = (1-s(e))*theta_start+s(e)*theta_goal;
p1 = p_robot(:,1);
p2 = p_robot(:,2);
p3 = stest(S1,theta(1))*[p_robot(:,3);1];
p4 = stest(S1,theta(1))*stest(S2,theta(2))*[p_robot(:,4);1];
p5 = stest(S1,theta(1))*stest(S2,theta(2))*stest(S3,theta(3))*[p_robot(:,5);1];
p6 = stest(S1,theta(1))*stest(S2,theta(2))*stest(S3,theta(3))*stest(S4,theta(4))*[p_robot(:,6);1];


p = [p1 p2 p3(1:3) p4(1:3) p5(1:3) p6(1:3)];
[~,nspheres]=size(p_robot);
combo = combnk(1:nspheres,2);
[ncombo,~]=size(combo);
[~,nobst]=size(p_obstacle);

for x = 1:ncombo
    
   d_s(x)=norm(p(:,combo(x,1))-p(:,combo(x,2)))-0.9-0.9; 
end

for y = 1:nspheres
    for z = 1:nobst
        d_o(y,z)=norm(p(:,y)-p_obstacle(:,z))-0.9-r_obstacle(z);
    end
end

if min(d_s) <= 0 || min(d_o(:))<=0
    c(e) = s(e);
else
    c(e)=0;
end

if 1 && all(c==0)
   A = 0;
else
    A=1;
end

end

end
