function [p, e, t] = triangulate(g, h, s)
  [p_b, e] = surface_mesh(g, h);
  n_b = size(p_b, 2);
  p_i = inner_nodes(g, h);
  p = [p_b, p_i];
  t = delaunay(p(1,:), p(2,:));
  t = orient(p, t);
  for k = 1:s
    p = smooth(p, t, n_b);
  end
  t = t';
