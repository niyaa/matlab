function [M_gamma, M, S] = fe_matrices2d(p, e, t, alpha, gamma, k)
  np = size(p, 2);
  nt = size(t, 2);
  switch (k)
   case 1
    dof_loc = 3;
    q = 2;
    shapefunction = @P1ShapeFunction;
    dof_glob = np;
   case 2
    dof_loc = 6;
    q = 4;
    shapefunction = @P2ShapeFunction;
    edges = GetEdges(np, e, t);
    TtoE = MakeTtoE(t, edges);
    dof_glob = np + size(edges, 2);
   otherwise
    disp('fe_matrices2d: approximation order not implemented!');
    return;
  end
  S = spalloc(dof_glob, dof_glob, 6*dof_glob); 
  M = spalloc(dof_glob, dof_glob, 6*dof_glob);
  M_gamma = spalloc(dof_glob, dof_glob, 6*dof_glob);
  
  for i = 1:nt 
    p1 = p(:,(t(1,i))); 
    p2 = p(:,(t(2,i)));
    p3 = p(:,(t(3,i)));
    [M_gamma_loc, M_loc, S_loc] = element_matrix2d(p1, p2, p3, alpha, ...
                             gamma, q, dof_loc, shapefunction);  

    if k == 1 
      idx = t(1:3,i); 
    else 
      idx = [t(1:3,i)' np+TtoE(1:3,i)']; 
    end
    S(idx, idx) = S(idx, idx) + S_loc; 
    M(idx, idx) = M(idx, idx) + M_loc;
    M_gamma(idx, idx) = M_gamma(idx, idx) + M_gamma_loc;
  end

