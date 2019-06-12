function A = min_norm(M,q_in0,q_in1)
    
    [~,npoints]=size(q_in1);
    
    for x = 1:npoints        
       dist(x) = norm(q_in0-M*q_in1(:,x)); 
    end
    
    [~,min_ind]=min(dist);
    
    A = min_ind;




end