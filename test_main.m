clear
clc

S = [1.00 0.00 -1.00 -1.00; 0.00 1.00 0.00 0.00; 0.00 0.00 0.00 0.00; 0.00 0.00 0.00 0.00; 0.00 0.00 0.00 0.00; 0.00 2.00 2.00 4.00];
M = [-1.00 0.00 0.00 -2.00; 0.00 0.00 1.00 4.00; 0.00 1.00 0.00 0.00; 0.00 0.00 0.00 1.00];
p_robot = [0.00 2.00 2.00 0.00 0.00 -2.00; 0.00 0.00 2.00 2.00 4.00 4.00; 0.00 0.00 0.00 0.00 0.00 0.00];
r_robot = [0.90 0.90 0.90 0.90 0.90 0.90];
p_obstacle = [2.57 0.98 2.24 0.38 4.39 2.15 2.04 1.83 0.63 4.46 -3.89 -1.52 2.98 5.09 5.20 5.68 0.41 -2.07 5.26 0.42 4.14 -2.62 3.72 -2.72 1.68 -5.78 3.40 3.24 1.26 4.45 -4.47 3.23 -1.86 5.75 0.95 3.33 0.05 -0.64 2.03 -0.31; 1.54 -0.74 -3.91 5.01 5.86 1.55 1.71 4.22 -2.74 -3.49 -3.47 4.02 -4.98 4.48 4.68 3.96 -2.01 0.86 -4.87 -5.26 -4.74 0.71 -2.63 4.35 -4.34 2.62 -1.20 -1.53 -1.01 1.16 -5.03 -2.61 1.33 0.42 2.34 3.21 5.24 -1.55 -5.69 -0.33; 2.65 4.27 -5.56 -5.44 5.51 -4.88 5.09 1.57 -0.30 2.28 -1.76 4.04 3.71 -2.94 -0.24 5.43 -5.36 0.94 -4.85 2.44 -5.86 -2.73 2.64 -3.77 -0.35 -4.63 -0.10 2.97 -5.51 2.78 0.86 3.49 -0.83 3.87 2.29 1.07 0.08 -0.82 -4.72 -4.11];
r_obstacle = [0.66 0.21 0.76 0.58 0.39 0.28 0.26 0.41 0.21 0.73 0.68 0.71 0.30 0.68 0.26 0.53 0.23 0.59 0.72 0.45 0.70 0.64 0.43 0.33 0.32 0.59 0.20 0.69 0.60 0.41 0.63 0.33 0.63 0.61 0.35 0.44 0.76 0.32 0.33 0.78];
theta_start = [-1.20; -1.52; 1.77; 0.13];
theta_goal = [0.10; -2.41; -2.86; 1.29];

[njoints,~]=size(theta_start);
forward_node.theta = theta_start;
forward_node.parent = 0;
T_forward{1}=forward_node;
backward_node.theta = theta_goal;
backward_node.parent = 0;
T_backward{1} = backward_node;


iter = 1;
while iter < 200
    %rand theta
    c_state = 1;
    while c_state==1
       rng shuffle
       theta_samp = 2*pi.*rand(njoints,1)-pi;
       c_state = coli_test(S,p_robot,r_robot,p_obstacle,r_obstacle,theta_start,theta_samp);
    end
    
    %f path
    [~,tfsize] = size(T_forward);
    for a = 1:tfsize
        
       d_f(a) = norm(theta_samp-T_forward{a}.theta); 
    end
    [~,minindf]=min(d_f);
    theta_fc = T_forward{minindf}.theta;
    theta_fc_ind = minindf;
    
    f_ccs = coli_test(S,p_robot,r_robot,p_obstacle,r_obstacle,theta_fc,theta_samp);
    
    if f_ccs == 0
       forward_node.theta = theta_samp;
       forward_node.parent = theta_fc_ind;
       T_forward = [T_forward forward_node];
        
    end
    
     %b path
    [~,tbsize] = size(T_backward);
    for b = 1:tbsize
        
       d_b(b) = norm(theta_samp-T_backward{b}.theta); 
    end
    [~,minindb]=min(d_b);
    theta_bc = T_backward{minindb}.theta;
    theta_bc_ind = minindb;
    
    b_ccs = coli_test(S,p_robot,r_robot,p_obstacle,r_obstacle,theta_bc,theta_samp);
    
    if b_ccs == 0
       backward_node.theta = theta_samp;
       backward_node.parent = theta_bc_ind;
       T_backward = [T_backward backward_node];
        
    end
    
    
    %sc
    iter=iter+1;
    [~,tbsize] = size(T_backward);
    [~,tfsize] = size(T_forward);
    
    for h = 1:tfsize
       f_cstate = isequal(ismember(theta_samp,T_forward{tfsize}.theta),ones(njoints,1)); 
    end
    
     for hh = 1:tbsize
       b_cstate = isequal(ismember(theta_samp,T_backward{tbsize}.theta),ones(njoints,1)); 
     end
    
    if f_cstate == b_cstate || f_cstate==1
        tfs = tfsize;
        tbs = tbsize;
        q=[];
        while tfs ~=0
           q = [T_forward{tfs}.theta q];
           fpar = T_forward{tfs}.parent;
           tfs = fpar;
        end
        if tbs ~= 1
            tbs = tbs-1;
        end
        while tbs ~=0
           q = [q T_backward{tbs}.theta];
           bpar = T_backward{tbs}.parent;
           tbs = bpar;
        end
        disp(mat2str(q,5))
        return
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end