function [ el2, eh1 ] = errors ( element_area, element_node, node_xy, ...
  f, element_num, element_order, nqe, node_num )

%*****************************************************************************80
%
%% ERRORS calculates the L2 and H1-seminorm errors.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    17 May 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real ELEMENT_AREA(ELEMENT_NUM), the area of elements.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
%    ELEMENT_NODE(I,J) is the global index of local node I in element J.
%
%    Input, real NODE_XY(2,NODE_NUM), the nodes.
%
%    Input, real F(NODE_NUM), the coefficients of the solution.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
%
%    Input, integer NQE, the number of points in the quadrature rule.
%    This is actually fixed at 13.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, real EL2, the L2 error.
%
%    Output, real EH1, the H1 seminorm error.
%
%  Local Parameters:
%
%    Local, real AR, the weight for a given quadrature point
%    in a given element.
%
%    Local, real BB, BX, BY, a basis function and its first
%    derivatives evaluated at a particular quadrature point.
%
%    Local, real EH1, the H1 seminorm error.
%
%    Local, real EL2, the L2 error.
%
%    Local, real UEX, UEXX, UEXY, the exact solution and its first
%    derivatives evaluated at a particular quadrature point.
%
%    Local, real UH, UHX, UHY, the computed solution and its first
%    derivatives evaluated at a particular quadrature point.
%
%    Local, real WQE(NQE), stores the quadrature weights.
%
%    Local, real X, Y, the coordinates of a particular
%    quadrature point.
%
%    Local, real XQE(NQE), YQE(NQE), stores the location
%    of quadrature points in a given element.
%
  el2 = 0.0E+00;
  eh1 = 0.0E+00;

  for element = 1 : element_num

    [ wqe, xqe, yqe ] = quad_e ( node_xy, element_node, element, ...
      element_num, element_order, node_num, nqe );

    for iq = 1 : nqe

      ar = element_area(element) * wqe(iq);
      x = xqe(iq);
      y = yqe(iq);

      uh = 0.0E+00;
      dudxh = 0.0E+00;
      dudyh = 0.0E+00;

      for in1 = 1 : element_order

        i = element_node(in1,element);

        [ bi, dbidx, dbidy ] = qbf ( x, y, element, in1, node_xy, ...
          element_node, element_num, element_order, node_num );

        uh    = uh    + bi    * f(i);
        dudxh = dudxh + dbidx * f(i);
        dudyh = dudyh + dbidy * f(i);

      end

      [ u, dudx, dudy ] = exact ( x, y );

      el2 = el2 + ( uh - u )^2 * ar;
      eh1 = eh1 + ( ( dudxh - dudx )^2 + ( dudyh - dudy )^2 ) * ar;

    end

  end

  el2 = sqrt ( el2 );
  eh1 = sqrt ( eh1 );

  fprintf ( 1, '\n' );
  fprintf ( 1, '*********************************************\n' );
  fprintf ( 1, '*                                           *\n' );
  fprintf ( 1, '*  ERRORS:                                  *\n' );
  fprintf ( 1, '*    L2 error =          %14f     *\n', el2 );
  fprintf ( 1, '*    H1-seminorm error = %14f     *\n', eh1 );
  fprintf ( 1, '*                                           *\n' );
  fprintf ( 1, '*********************************************\n' );

  return
end