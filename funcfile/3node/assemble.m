function [ a, f ] = assemble ( node_num, node_xy, element_order, ...
  element_num, element_node, nq, wq, xq, yq, element_area )


  if ( element_order == 3 )
    nz_max = 7 * node_num;
  else
    nz_max = 19 * node_num;
  end

   %a = sparse ( [], [], [], node_num, node_num, nz_max );
   a=zeros(node_num);
%
%  Initialize the arrays to zero.
%
  f(1:node_num) = 0.0;
%
%  The actual values of A and F are determined by summing up
%  contributions from all the elements.
%
  for element = 1 : element_num
%
%  Consider a quadrature point QUAD, with coordinates (X,Y).
%
    for quad = 1 : nq

      x = xq(quad,element);
      y = yq(quad,element);
      w = element_area(element) * wq(quad);
%
%  Consider one of the basis functions, which will play the
%  role of test function in the integral.
%
      for test = 1 : element_order

        i = element_node(test,element);

        [ bi, dbidx, dbidy ] = qbf ( x, y, element, test, node_xy, ...
          element_node, element_num, element_order, node_num );

        f(i) = f(i) + w * rhs ( x, y ) * bi;
%
%  Consider a basis function used to form the value of the solution function.
%
%  Logically, this term goes in entry A(I,J).  Because of the
%  band matrix storage, entry (I,J) is actually stored in
%  A(I-J+2*NHBA+1,J).
%
        for basis = 1 : element_order

          j = element_node(basis,element);

          [ bj, dbjdx, dbjdy ] = qbf ( x, y, element, basis, node_xy, ...
            element_node, element_num, element_order, node_num );

          aij = dbidx * dbjdx + dbidy * dbjdy;

          a(i,j) = a(i,j) + w * aij;

        end

      end

    end

  end

  return
end