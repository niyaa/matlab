function [p_b, e] = surface_mesh(g, h)
  n_bs = size(g,2); % Anzahl Randsegmente
  g(:, n_bs+1) = g(:, 1);
  x = [];
  y = [];
  e3 = [];
  n = 0;
  for k = 1:n_bs
    H = norm([g(1,k+1) - g(1,k); g(2,k+1) - g(2,k)]);
    n_el = ceil(H/h);
    nodes_x = linspace(g(1,k), g(1,k+1), n_el+1);
    nodes_y = linspace(g(2,k), g(2,k+1), n_el+1);
    x = [x nodes_x(1:end-1)]; 
    y = [y nodes_y(1:end-1)];
    e3 = [e3 k*ones(1, n_el)];
  end
  n = length(x);
  p_b = zeros(2, n);
  p_b(1, :) = x;
  p_b(2, :) = y;
  e = [1:n; [2:n 1]; e3];
  
  
  
  
  