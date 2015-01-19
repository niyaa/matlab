function TtoE = MakeTtoE(t, edges)
% Creates a triangle to edge connectivity matrix.
%
% INPUT:   t,edges  mesh data
%
% OUTPUT:  TtoE  3 by nt matrix, nt the number of elements. The i-th
%                column contains the indices (with respect to edges) of the
%                edges of the triangle i (with respect to t). 

  nt = size(t, 2);
  TtoE = zeros(3, nt);
  e1 = edges(1,:);
  e2 = edges(2,:);
  j2vec = [2 3 1];
  for k = 1:nt
    for j1 = 1:3
      j2 = j2vec(j1);
      tj1k = t(j1,k);
      tj2k = t(j2,k);
      TtoE(j1,k) = find((e1 == tj1k | e1 == tj2k) & (e2 == tj1k | e2 == tj2k));
    end
  end
	
	