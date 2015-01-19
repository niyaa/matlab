function u = fe_solution2d(p,e,t,alpha,gamma,f,bc,g,k)

[M,S,b] = fe_matrices2d(p,e,t,alpha,gamma,f,k);
A = M + S;
ne = size(e,2);
np = size(p,2);

switch k
    case 1
        for i = 1:ne
            if bc(e(3,i)) == 'n'                              
                b_loc = fe_surfaceint(p(:,e(1,i)),p(:,e(2,i)),g,e(3,i),k);
                b(e(1,i)) = b(e(1,i)) + b_loc(1);
                b(e(2,i)) = b(e(2,i)) + b_loc(2);
            end    
        end
        for i = 1:ne
            if bc(e(3,i)) == 'd'                              
                b(e(1,i)) = feval(g,p(1,e(1,i)),p(2,e(1,i)),e(3,i));
                b(e(2,i)) = feval(g,p(1,e(2,i)),p(2,e(2,i)),e(3,i));
                A(e(1,i),:) = 0;
                A(e(2,i),:) = 0;
                A(e(1,i),e(1,i)) = 1;
                A(e(2,i),e(2,i)) = 1;
            end    
        end          
              
    case 2
        for i = 1:ne
            if bc(e(3,i)) == 'n'                              
                b_loc = fe_surfaceint(p(:,e(1,i)),p(:,e(2,i)),g,e(3,i),k);
                b(e(1,i)) = b(e(1,i)) + b_loc(1);
                b(e(2,i)) = b(e(2,i)) + b_loc(2);
                b(np+i) = b(np+i) + b_loc(3);
            end    
        end
        for i = 1:ne
            if bc(e(3,i)) == 'd'                              
                p_med = 0.5*(p(:,e(1,i))+p(:,e(2,i))); 
                b(e(1,i)) = feval(g,p(1,e(1,i)),p(2,e(1,i)),e(3,i));
                b(e(2,i)) = feval(g,p(1,e(2,i)),p(2,e(2,i)),e(3,i));
                b(np+i) = feval(g,p_med(1),p_med(2),e(3,i));
                A(e(1,i),:) = 0;
                A(e(2,i),:) = 0;
                A(np+i,:) = 0;
                A(e(1,i),e(1,i)) = 1;
                A(e(2,i),e(2,i)) = 1;
                A(np+i,np+i) = 1;
            end    
        end
        
end

u = A\b;