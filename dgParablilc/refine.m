function [p, e, t] = refine(p, e, t)
  np = size(p, 2);
  nt = size(t, 2);
  ne = size(e, 2);
  
  edges = GetEdges(np, e, t);
  p1_idx = edges(1, :);
  p2_idx = edges(2, :);
  pnew = 0.5*(p(:, p1_idx) + p(:, p2_idx));
  p = [p pnew];
  
  enew = zeros(3, 2*ne);
  enew(1, 1:ne) = e(1, :);
  enew(1, ne+1:end) = np+1:np+ne;
  enew(2, ne+1:end) = e(2, :);
  enew(2, 1:ne) = np+1:np+ne;
  enew(3, 1:ne) = e(3, :);
  enew(3, ne+1:end) = e(3, :);
  e = enew;
  
  TtoE = MakeTtoE(t, edges);
  tnew1 = [];
  tnew2 = [];
  tnew3 = [];
  for k = 1:nt
    father = t(1:3, k);
    newnodes = np + TtoE(1:3, k);
    tnew1 = [tnew1 newnodes(1) father(1)   father(2)   father(3)  ]; 
    tnew2 = [tnew2 newnodes(2) newnodes(1) newnodes(2) newnodes(3)];
    tnew3 = [tnew3 newnodes(3) newnodes(3) newnodes(1) newnodes(2)];
  end  
  t = [tnew1; tnew2; tnew3];
  
  
  
  