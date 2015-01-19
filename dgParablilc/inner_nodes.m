function p_i = inner_nodes(g, h)
  n_bs = size(g,2); % Anzahl Randsegmente
  x = [];
  y = [];
  llx = min(g(1,:));
  lly = min(g(2,:));
  urx = max(g(1,:));
  ury = max(g(2,:));
  n_elx = ceil((urx - llx)/h);
  n_ely = ceil((ury - lly)/h);
  nodes_x = linspace(llx, urx, n_elx+1);
  nodes_y = linspace(lly, ury, n_ely+1);
  for k = 2:n_ely
    x = [x nodes_x(2:n_elx)];
    y = [y nodes_y(k)*ones(1,n_elx-1)];
  end
  in = my_inpolygon(x, y, g(1,:), g(2,:));
  x = x(in);
  y = y(in);
  n = length(x);
  p_i = zeros(2, n);
  if (~isempty(p_i))
    p_i(1, :) = x;
    p_i(2, :) = y;
  end
  
  
  
  
  
  