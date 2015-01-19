function [ a, f ] = assemble ( node_num, node_xy, element_order, ...
  element_num, element_node, nq, wq, xq, yq, element_area )

%*****************************************************************************80
%
%% ASSEMBLE assembles the coefficient matrix A and right hand side F.
%
%  Discussion:
%
%    The matrix is known to be banded.  A special matrix storage format
%    is used to reduce the space required.  Details of this format are
%    discussed in the routine DGB_FA.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 March 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, real NODE_XY(2,NODE_NUM), the nodes.
%
%    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
%    ELEMENT_NODE(I,J) is the global index of local node I in element J.
%
%    Input, integer NQ, the number of quadrature points used in assembly.
%
%    Input, real WQ(NQ), quadrature weights.
%
%    Input, real XQ(NQ,ELEMENT_NUM), YQ(NQ,ELEMENT_NUM), the
%    coordinates of the quadrature points in each element.
%
%    Input, real ELEMENT_AREA(ELEMENT_NUM), the area of elements.
%
%    Input, integer IB, the half-bandwidth of the matrix.
%
%    Output, real A(3*IB+1,NODE_NUM), the NODE_NUM by NODE_NUM
%    coefficient matrix, stored in a compressed format.
%
%    Output, real F(NODE_NUM), the right hand side.
%
%  Local parameters:
%
%    Local, real BB, BX, BY, the value of some basis function
%    and its first derivatives at a quadrature point.
%
%    Local, real BBB, BBX, BBY, the value of another basis
%    function and its first derivatives at a quadrature point.
%

%
%  Use MATLAB's sparse matrix storage format.
%
  if ( element_order == 3 )
    nz_max = 7 * node_num;
  else
    nz_max = 19 * node_num;
  end

  a = sparse ( [], [], [], node_num, node_num, nz_max );
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