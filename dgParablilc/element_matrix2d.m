function [M_gamma, M, S] = element_matrix2d(p1, p2, p3, alpha, gamma, ...
                          q, dof, shapefunction)

  [loc, glob, weights, DFinv] = quad_nodes(p1, p2, p3, q);
  M_gamma = zeros(dof); 
  M = zeros(dof); 
  S = zeros(dof); 
  for k = 1:length(weights)
    x_k = glob(1, k);
    y_k = glob(2, k);
    w_k = weights(k);
    A_k = alpha(x_k, y_k);
    c_k = gamma(x_k, y_k);
    for i = 1:dof
      [v_i, sgrad_i] = shapefunction(i, loc(:, k));
      grad_i = DFinv*sgrad_i;
      Agrad_i = A_k*grad_i;
      cv_i = c_k*v_i;
      for j = i:dof
        [v_j, sgrad_j] = shapefunction(j, loc(:, k));
        grad_j = DFinv*sgrad_j;
        M_gamma(i, j) = M_gamma(i, j) + w_k*cv_i*v_j;
        M(i, j) = M(i, j) + w_k*v_i*v_j;
        S(i, j) = S(i, j) + w_k*grad_j'*Agrad_i;
        M_gamma(j,i) = M_gamma(i,j);
        M(j,i) = M(i,j);
        S(j,i) = S(i,j);
      end
    end
  end
 
