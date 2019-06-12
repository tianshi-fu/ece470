clear
clc

%path planning

S = [0.00 1.00; 1.00 0.00; 0.00 0.00; 2.00 0.00; 0.00 -2.00; 0.00 -2.00];
M = [0.00 0.00 1.00 0.00; 0.00 -1.00 0.00 2.00; 1.00 0.00 0.00 0.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 0.00 0.00 0.00; 0.00 0.00 2.00 2.00; 0.00 -2.00 -2.00 0.00];
r_robot = [0.90 0.90 0.90 0.90];
p_obstacle = [4.58 2.64 4.44 0.47 3.93 3.20 -2.02 -3.17 3.47 2.56 -2.00 -0.17 1.30 -0.24 4.05 -0.10 -0.65 2.35 -2.25 -4.08 2.35 -3.87 3.08 -4.89 3.30 -0.99 4.00; -0.04 -1.61 2.78 -2.56 4.00 3.34 -2.24 -1.36 4.76 -1.46 -0.72 0.79 -2.22 -4.27 -3.94 -1.08 4.61 2.06 -4.11 -1.83 2.58 0.91 -4.20 3.24 -2.45 -3.97 4.97; -2.13 2.00 -3.96 4.90 1.43 3.39 3.98 -4.65 -2.70 4.86 4.92 -4.50 -4.41 -2.64 0.09 2.23 4.72 4.94 -4.82 0.03 -3.27 3.78 -3.58 3.56 -3.76 1.46 0.87];
r_obstacle = [3.28 1.47 3.17 2.58 1.78 1.56 1.28 2.44 2.12 1.37 2.70 1.21 1.88 0.69 4.44 0.90 3.26 4.79 3.39 0.80 1.17 3.52 1.00 3.03 2.03 2.93 3.72];
theta_start = [2.01; -1.56];
theta_goal = [-0.14; -0.83];

forward_node.theta = theta_start;
forward_node.parent = 0;
T_forward{1} = forward_node;
backward_node.theta = theta_goal;
backward_node.parent = 0;
T_backward{1} = backward_node;


iter = 1;
while iter < 200
    %generate random theta
    c_state = 1;
    while c_state ==1 %loop until a theta in free space is found
        rng shuffle
        theta_samp = (pi-(-pi)).*rand(2,1) + (-pi); %create random theta array from -pi to pi
        c_state = colli_2j(S,theta_start,theta_samp,p_robot,r_robot,p_obstacle,r_obstacle); %check if the random theta is in free space
    end
    
    %forward path
    [~,Tf_size]=size(T_forward);
    for k = 1:Tf_size
    dist_f(k) = norm(theta_samp-T_forward{k}.theta); %find distances between theta_samp and all thetas stored in forward tree
    end
    [~,min_indf] = min(dist_f);
    theta_fclosest = T_forward{min_indf}.theta; %find cloeset theta in forward tree to theta_samp
    theta_fclosest_ind = min_indf; %find the index of cloeset theta in forward tree to theta_samp 
    
    forward_cstate = colli_2j(S,theta_fclosest,theta_samp,p_robot,r_robot,p_obstacle,r_obstacle); %check if the path from theta_fclosest to theta_samp is collisioni free
    
    if forward_cstate == 0
        forward_node.theta = theta_samp; %assign theta_samp to .theta data structure
        forward_node.parent = theta_fclosest_ind; %assign theta_fclosest_ind to .parent data structure 
        T_forward = [T_forward forward_node]; %add the new node to forward tree
    end
    
    %backward path
    [~,Tb_size]=size(T_backward);
    for m = 1:Tb_size
    dist_b(m) = norm(theta_samp-T_backward{m}.theta); %find distances between theta_samp and all thetas stored in backward tree
    end
    [~,min_indb] = min(dist_b);
    theta_bclosest = T_backward{min_indb}.theta; %find cloeset theta in backward tree to theta_samp
    theta_bclosest_ind = min_indb; %find the index of cloeset theta in backward tree to theta_samp 
    
    backward_cstate = colli_2j(S,theta_bclosest,theta_samp,p_robot,r_robot,p_obstacle,r_obstacle); %check if the path from theta_fclosest to theta_samp is collisioni free
    
    if backward_cstate == 0
        backward_node.theta = theta_samp; %assign theta_samp to .theta data structure
        backward_node.parent = theta_bclosest_ind; %assign theta_bclosest_ind to .parent data structure 
        T_backward = [T_backward backward_node]; %add the new node to backward tree
    end
    
    %terminating condition
     iter = iter+1;
     [~,Tf_size]=size(T_forward);
     [~,Tb_size]=size(T_backward);
     
     for w = 1:Tf_size
     f_state = isequal(ismember(theta_samp,T_forward{w}.theta),[1;1]);
     end
     for ww = 1:Tb_size
     b_state = isequal(ismember(theta_samp,T_backward{ww}.theta),[1;1]);
     end
     
     if f_state == b_state
         q = [];
         while Tf_size ~= 0     
            q = [T_forward{Tf_size}.theta q];
            fpar = T_forward{Tf_size}.parent;
            Tf_size = fpar;
         end    
         Tb_size = Tb_size -1;
         while Tb_size ~=0
             q = [q T_backward{Tb_size}.theta ];
             bpar = T_backward{Tb_size}.parent;
             Tb_size = bpar;
         end
         disp(mat2str(q,5))
            return
     end
end
