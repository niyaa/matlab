function edges = GetEdges(np, e, t)
% Locates the edges of a finite element mesh.
%
% INPUT:   np  number of np vertices
%          e   boundary edges
%          t   elements
%
% OUTPUT:  edges  2 by ne matrix, ne the number of edges. The i-th
%                 column contains the indices (with respect to p) of the
%                 edge endpoints
%
% WARNING These are essentially lines taken from the refinemesh
%         routine of the PDE toolbox. Be sure to own a license!!!
  nt = size(t,2);
  it = 1:nt;
  ip1 = t(1,it);
  ip2 = t(2,it);
  ip3 = t(3,it);
  A = sparse(ip1, ip2, -1, np, np);
  A = A + sparse(ip2, ip3, -1, np, np);
  A = A + sparse(ip3, ip1, -1, np, np);
  A = -((A + A.') < 0);
  ie = full(A(e(1,:) + (e(2,:)-1)*np)) == -1;
  ie = find(ie);                            
  ip = (np + 1):(np + length(ie));
  A(e(1,ie) + np*(e(2,ie)-1)) = ip;
  A(e(2,ie) + np*(e(1,ie)-1)) = ip;
  [i1,i2] = find(A == -1 & A.' == -1);
  i = find(i2 > i1);
  i1 = i1(i);
  i2 = i2(i);
  edges = zeros(2, length(ie) + length(i1));
  edges(1, 1:length(ie)) = e(1,ie);
  edges(2, 1:length(ie)) = e(2,ie);
  edges(1, length(ie)+1:size(edges,2)) = i1';
  edges(2, length(ie)+1:size(edges,2)) = i2';

