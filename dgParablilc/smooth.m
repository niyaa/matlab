function p = smooth(p, tri, n_b);
  n_p = size(p, 2);
  x = p(1,:); y = p(2,:);
  for k = n_b+1:n_p
    [idx, help] = find(tri == k);
    n_l = length(idx);
    x_new = 0;
    y_new = 0;
    for j = 1:n_l;
      local = mod(find(tri(idx(j), :) == k), 3) + 1;
      x_new = x_new + x(tri(idx(j), local));
      y_new = y_new + y(tri(idx(j), local));
    end
    x(k) = x_new/n_l;
    y(k) = y_new/n_l;
  end
  p = [x; y];
  