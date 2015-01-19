function [P,rhs] = apply_bc(P_old, rhs_old, p, e, t, bc, g, k, time_deg, time, u0);

% Applies boundary conditions to the system matrix and the right hand side

ne = size(e,2);
np = size(p,2);

switch time_deg
    case 0
        P = P_old;
        rhs = rhs_old;
                
        switch k
            case 1
                for i = 1:ne
                    if bc(e(3,i)) == 'n'                              
                        rhs_loc = fe_surfaceint(p(:,e(1,i)), p(:,e(2,i)), g, e(3,i), k, time_deg, time);
                        rhs(e(1:2,i)) = rhs(e(1:2,i)) + rhs_loc;                     
                    end    
                end
                for i = 1:ne
                    if bc(e(3,i)) == 'd'                              
                        rhs(e(1,i)) = g(p(1,e(1,i)), p(2,e(1,i)), time(2), e(3,i));
                        rhs(e(2,i)) = g(p(1,e(2,i)), p(2,e(2,i)), time(2), e(3,i));
                        P(e(1:2,i),:) = 0;
                        P(e(1,i),e(1,i)) = 1;
                        P(e(2,i),e(2,i)) = 1;
                    end    
                end          
              
            case 2
                for i = 1:ne
                    if bc(e(3,i)) == 'n'                              
                        rhs_loc = fe_surfaceint(p(:,e(1,i)), p(:,e(2,i)), g, e(3,i), k, time_deg, time);
                        rhs(e(1:2,i)) = rhs(e(1:2,i)) + rhs_loc(1:2);
                        rhs(np+i) = rhs(np+i) + rhs_loc(3);
                    end    
                end
                for i = 1:ne
                    if bc(e(3,i)) == 'd'                              
                        p_mean = 0.5*(p(:,e(1,i))+p(:,e(2,i)));
                        rhs(e(1,i)) = g(p(1,e(1,i)), p(2,e(1,i)), time(2), e(3,i));
                        rhs(e(2,i)) = g(p(1,e(2,i)), p(2,e(2,i)), time(2), e(3,i));
                        rhs(np+i) = g(p_mean(1), p_mean(2), time(2), e(3,i));
                        P(e(1:2,i),:) = 0;
                        P(np+i,:) = 0;
                        P(e(1,i),e(1,i)) = 1;
                        P(e(2,i),e(2,i)) = 1;
                        P(np+i,np+i) = 1;
                    end    
                end       
        end
        
    case 1
        P = P_old;
        block_size = 0.5*size(P,1);
        rhs = rhs_old;
        
        switch k
            case 1
                for i = 1:ne
                    if bc(e(3,i)) == 'n'                              
                        rhs_loc = fe_surfaceint(p(:,e(1,i)), p(:,e(2,i)), g, e(3,i), k, time_deg, time);
                        
                        rhs(e(1:2,i)) = rhs(e(1:2,i)) + rhs_loc(1:2);                     % block F1
                        
                        rhs(block_size + e(1:2,i)) = rhs(block_size + e(1:2,i)) + rhs_loc(1:2); % block F2
                    end    
                end
                for i = 1:ne
                    if bc(e(3,i)) == 'd'                              
                        rhs(e(1,i)) = u0(e(1,i));  % block F1
                        rhs(e(2,i)) = u0(e(2,i));
                        
                        rhs(block_size + e(1,i)) = g(p(1,e(1,i)), p(2,e(1,i)), time(2), e(3,i))-u0(e(1,i)); % block F2   
                        rhs(block_size + e(2,i)) = g(p(1,e(2,i)), p(2,e(2,i)), time(2), e(3,i))-u0(e(2,i));              
                        
                        P(e(1:2,i),:) = 0;                                           % block P1 Q1
                        P(e(1,i),e(1,i)) = 1;
                        P(e(2,i),e(2,i)) = 1;
                        
                        P(block_size + e(1:2,i),:) = 0;                              % block P2 Q2
                        P(block_size + e(1,i), block_size + e(1,i)) = 1;
                        P(block_size + e(2,i), block_size + e(2,i)) = 1;
                        
                  end    
                end          
              
            case 2
                
                for i = 1:ne
                    if bc(e(3,i)) == 'n'
                        rhs_loc = fe_surfaceint(p(:,e(1,i)), p(:,e(2,i)), g, e(3,i), k, time_deg, time);
                        
                        rhs(e(1:2,i)) = rhs(e(1:2,i)) + rhs_loc(1:2);               % block F1
                        rhs(np+i) = rhs(np+i) + rhs_loc(3);

                        rhs(block_size + e(1:2,i)) = rhs(block_size + e(1:2,i)) + rhs_loc(1:2);  % block F2
                        rhs(block_size + np+i) = rhs(block_size + np+i) + rhs_loc(3);
                    end    
                end
                for i = 1:ne
                    if bc(e(3,i)) == 'd'                              
                        p_mean = 0.5*(p(:,e(1,i))+p(:,e(2,i)));
                         
                        rhs(e(1,i)) = u0(e(1,i));  % block F1
                        rhs(e(2,i)) = u0(e(2,i));
                        rhs(np+i) = u0(np + i);
			
                        rhs(block_size + e(1,i)) = g(p(1,e(1,i)), p(2,e(1,i)), time(2), e(3,i))-u0(e(1,i)); % block F2   
                        rhs(block_size + e(2,i)) = g(p(1,e(2,i)), p(2,e(2,i)), time(2), e(3,i))-u0(e(2,i));              
                        rhs(block_size + np+i) = g(p_mean(1), p_mean(2), time(2), e(3,i)) - u0(np + i);
                        
                        P(e(1:2,i),:) = 0;                                           % block P1 Q1
                        P(np+i,:) = 0;
                        P(e(1,i),e(1,i)) = 1;
                        P(e(2,i),e(2,i)) = 1;
                        P(np+i,np+i) = 1;
                        
                        P(block_size + e(1:2,i),:) = 0;                              % block P2 Q2
                        P(block_size + np+i,:) = 0;
                        P(block_size + e(1,i), block_size + e(1,i)) = 1;
                        P(block_size + e(2,i), block_size + e(2,i)) = 1;
                        P(block_size + np+i, block_size + np+i) = 1;
                    end    
                end   
        end
end