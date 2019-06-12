clear
clc

%path planning

S = [0.00 0.00 0.00 -1.00 0.00 0.00 0.00; 0.00 0.00 -1.00 0.00 0.00 0.00 0.00; 0.00 0.00 0.00 0.00 1.00 0.00 1.00; 0.00 -1.00 -4.00 0.00 0.00 0.00 0.00; -1.00 0.00 0.00 4.00 4.00 -1.00 2.00; 0.00 0.00 2.00 0.00 0.00 0.00 0.00];
M = [0.00 -1.00 0.00 -4.00; 0.00 0.00 1.00 0.00; -1.00 0.00 0.00 0.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 0.00 0.00 -2.00 -4.00 -4.00 -2.00 -2.00 -4.00; 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00; 0.00 -2.00 -4.00 -4.00 -4.00 -2.00 -2.00 0.00 0.00];
r_robot = [0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90];
p_obstacle = [1.44 -2.72 4.18 0.51 4.04 1.07 4.79 4.44 -4.03 3.36 1.83 -0.28 -1.50 3.06 3.38; 4.21 4.92 4.42 3.16 -3.91 4.83 0.67 -1.67 2.34 -0.17 -4.89 1.86 4.49 1.55 -1.51; 1.40 -1.99 2.48 -1.95 2.51 -4.07 -3.55 4.62 3.18 3.95 4.43 1.46 -4.22 2.31 3.86];
r_obstacle = [3.65 0.82 1.78 0.76 4.31 1.64 2.58 0.57 3.69 1.11 2.72 0.76 0.98 2.38 2.33];
theta_start = [0.37; 0.33; -0.59; -2.52; 1.72; -2.11; 1.82];
theta_goal = [-2.26; 1.53; -1.64; 3.05; 1.33; -2.20; 1.61];

forward_node.theta = theta_start;
forward_node.parent = 0;
T_forward{1} = forward_node;
backward_node.theta = theta_goal;
backward_node.parent = 0;
T_backward{1} = backward_node;
[njoints,~] = size(theta_start);

iter = 1;
while iter < 200
    %generate random theta
    c_state = 1;
    while c_state ==1 %loop until a theta in free space is found
        rng shuffle
        theta_samp = (pi-(-pi)).*rand(njoints,1) + (-pi); %create random theta array from -pi to pi
        c_state = colli_nj(S,theta_start,theta_samp,p_robot,r_robot,p_obstacle,r_obstacle); %check if the random theta is in free space
    end
    
    %forward path
    [~,Tf_size]=size(T_forward);
    for k = 1:Tf_size
    dist_f(k) = norm(theta_samp-T_forward{k}.theta); %find distances between theta_samp and all thetas stored in forward tree
    end
    [~,min_indf] = min(dist_f);
    theta_fclosest = T_forward{min_indf}.theta; %find cloeset theta in forward tree to theta_samp
    theta_fclosest_ind = min_indf; %find the index of cloeset theta in forward tree to theta_samp 
    
    forward_cstate = colli_nj(S,theta_fclosest,theta_samp,p_robot,r_robot,p_obstacle,r_obstacle); %check if the path from theta_fclosest to theta_samp is collisioni free
    
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
    
    backward_cstate = colli_nj(S,theta_bclosest,theta_samp,p_robot,r_robot,p_obstacle,r_obstacle); %check if the path from theta_fclosest to theta_samp is collisioni free
    
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
     f_state = isequal(ismember(theta_samp,T_forward{w}.theta),ones(njoints,1));
     end
     for ww = 1:Tb_size
     b_state = isequal(ismember(theta_samp,T_backward{ww}.theta),ones(njoints,1));
     end
     
     if f_state == b_state && f_state == 1
         q = [];
         tfs = Tf_size;
         tbs = Tb_size;
         while tfs ~= 0     
            q = [T_forward{tfs}.theta q];
            fpar = T_forward{tfs}.parent;
            tfs = fpar;
         end    
         if tbs ~= 1
         tbs = tbs -1;
         end
         while tbs ~=0
             q = [q T_backward{tbs}.theta];
             bpar = T_backward{tbs}.parent;
             tbs = bpar;
         end
         disp(mat2str(q,5))
            return
     end
     disp(iter)
end
