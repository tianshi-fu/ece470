function [A,B,C]=calc_b_j(M,q_in0,q_in1)

      ident = eye(6);

      b = q_in0 - M*q_in1; 
      C = norm(q_in0 - M*q_in1)^2;
      
      
      
      J1 = ss61(ident(:,1))*M*q_in1;
      J2 = ss61(ident(:,2))*M*q_in1;
      J3 = ss61(ident(:,3))*M*q_in1;
      J4 = ss61(ident(:,4))*M*q_in1;
      J5 = ss61(ident(:,5))*M*q_in1;
      J6 = ss61(ident(:,6))*M*q_in1;
       
      J = [J1 J2 J3 J4 J5 J6];
   
   
      A = J' * J;
      B = J' * b;
   


end