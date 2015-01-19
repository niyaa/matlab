function tri = orient(p, tri)
%ORIENT Ensures consistent orientation of a triangulation.
%   Usage:  TRI = ORIENT(P, TRI) 
%   Input:  P    (2 by n_p) coordinate matrix
%           TRI  (n_t by 3) element index matrix
%   
%   Output: TRI  consistently oriented element index matrix
%
%   The routine checks the orientation of each element 
%   by calculating the determinant of the mapping from the 
%   reference element. In case of a negative orientation, 
%   the index ordering is reversed. 
  x = p(1,:);
  y = p(2,:);
  x21 = x(tri(:, 2)) - x(tri(:, 1));
  x31 = x(tri(:, 3)) - x(tri(:, 1));
  y21 = y(tri(:, 2)) - y(tri(:, 1));
  y31 = y(tri(:, 3)) - y(tri(:, 1));
  det = x21.*y31 - x31.*y21;
  neg = find(det < 0);
  dum = tri(neg, 1);
  tri(neg, 1) = tri(neg, 2);
  tri(neg, 2) = dum;
  
  
  
  
  
