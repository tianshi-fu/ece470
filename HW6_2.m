clear
clc
%path planning rrt?

S = [0.00 0.00; 1.00 0.00; 0.00 0.00; 0.00 0.00; 0.00 0.00; 2.00 -1.00];
M = [0.00 0.00 1.00 2.00; 1.00 0.00 0.00 4.00; 0.00 1.00 0.00 0.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 2.00 2.00 2.00; 0.00 0.00 2.00 4.00; 0.00 0.00 0.00 0.00];
r_robot = [0.90 0.90 0.90 0.90];
p_obstacle = [-2.91 1.06 3.57 -1.55 4.44 -1.73 3.82 -4.89 2.39 -3.74 2.59 -3.65 -0.99 -3.63 0.74; -3.05 -2.70 -3.42 4.91 4.50 -4.64 -0.88 0.73 -4.33 3.79 2.12 2.47 4.48 4.48 1.48; -2.75 -1.97 1.08 4.77 -2.94 -3.14 2.62 4.80 2.88 4.04 3.24 -0.53 -4.74 -3.48 4.13];
r_obstacle = [1.61 1.86 1.01 1.47 1.81 3.18 2.16 2.56 1.93 3.51 0.50 2.80 0.81 0.64 0.53];
theta_start = [-1.16; 0.30];
theta_goal = [0.38; 2.83];


S1 = S(:,1);
S2 = S(:,2);

max_iter = 1000;
i = 1;
done = 0;

q_start = theta_start;
Tree = q_start;

while i <= max_iter
    
    q_samp = (pi-(-pi)).*rand(2,1) + (-pi); %sample random point that has limit from -pi to pi
    
    [~,T_size]=size(Tree); 
    for a = 1:T_size %find nearest point in T to q_samp
    q_nearest = min(norm(q_samp-Tree(:,a)));
    end
    
    D = norm(q_samp-q_nearest); %distance between q_samp and q_nearest
    esp = 0.1; %resolution
    n = ceil(D/esp)+1; %n samples
    s = linspace(0,1,n); %create fractional distances
    
    for b = 1:n
        theta = (1-s(b))*q_nearest + s(b) * q_samp;
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
            c_status = 1;
            break   %breaks out of path collision checking
        else
            c_status = 0;
        end
    end
    if c_status == 0
    q_new = q_samp;
           Tree = [Tree, q_new];  %add q_samp to T
           
           if isequal(q_new,theta_goal) == 1
               done = 1;
               disp('done!')
               return
           end
    end
    disp(i) 
    i = i + 1;
end















