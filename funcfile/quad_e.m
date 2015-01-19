function [ wqe, xqe, yqe ] =  quad_e ( node_xy, element_node, ...
  element, element_num, element_order, node_num, nqe )

%*****************************************************************************80
%
%% QUAD_E sets up quadrature information for a 13-point rule in a given element.
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
%    Input, real NODE_XY(2,NODE_NUM), the nodes.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
%    ELEMENT_NODE(I,J) is the global index of local node I in element J.
%
%    Input, integer ELEMENT, the index of the element for which the quadrature
%    points are to be computed.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer NQE, the number of points in the quadrature rule.
%    This is actually fixed at 13.
%
%    Output, real WQE(13), the quadrature weights.
%
%    Output, real XQE(NQE), YQE(NQE), the X and Y coordinates
%    of the quadrature points.
%
  wqe = zeros ( 13, 1 );

  for i = 1 : 3
    wqe(i) = 0.175615257433204;
    ii = i + 3;
    wqe(ii) = 0.053347235608839;
    ii = i + 6;
    iii = ii + 3;
    wqe(ii) = 0.077113760890257;
    wqe(iii) = wqe(ii);
  end

  wqe(13) = -0.14957004446767;

  z1 = 0.479308067841923;
  z2 = 0.260345966079038;
  z3 = 0.869739794195568;
  z4 = 0.065130102902216;
  z5 = 0.638444188569809;
  z6 = 0.312865496004875;
  z7 = 0.048690315425316;

  ip1 = element_node(1,element);
  ip2 = element_node(2,element);
  ip3 = element_node(3,element);
  x1 = node_xy(1,ip1);
  x2 = node_xy(1,ip2);
  x3 = node_xy(1,ip3);
  y1 = node_xy(2,ip1);
  y2 = node_xy(2,ip2);
  y3 = node_xy(2,ip3);

  xqe( 1) = z1 * x1 + z2 * x2 + z2 * x3;
  yqe( 1) = z1 * y1 + z2 * y2 + z2 * y3;
  xqe( 2) = z2 * x1 + z1 * x2 + z2 * x3;
  yqe( 2) = z2 * y1 + z1 * y2 + z2 * y3;
  xqe( 3) = z2 * x1 + z2 * x2 + z1 * x3;
  yqe( 3) = z2 * y1 + z2 * y2 + z1 * y3;
  xqe( 4) = z3 * x1 + z4 * x2 + z4 * x3;
  yqe( 4) = z3 * y1 + z4 * y2 + z4 * y3;
  xqe( 5) = z4 * x1 + z3 * x2 + z4 * x3;
  yqe( 5) = z4 * y1 + z3 * y2 + z4 * y3;
  xqe( 6) = z4 * x1 + z4 * x2 + z3 * x3;
  yqe( 6) = z4 * y1 + z4 * y2 + z3 * y3;
  xqe( 7) = z5 * x1 + z6 * x2 + z7 * x3;
  yqe( 7) = z5 * y1 + z6 * y2 + z7 * y3;
  xqe( 8) = z5 * x1 + z7 * x2 + z6 * x3;
  yqe( 8) = z5 * y1 + z7 * y2 + z6 * y3;
  xqe( 9) = z6 * x1 + z5 * x2 + z7 * x3;
  yqe( 9) = z6 * y1 + z5 * y2 + z7 * y3;
  xqe(10) = z6 * x1 + z7 * x2 + z5 * x3;
  yqe(10) = z6 * y1 + z7 * y2 + z5 * y3;
  xqe(11) = z7 * x1 + z5 * x2 + z6 * x3;
  yqe(11) = z7 * y1 + z5 * y2 + z6 * y3;
  xqe(12) = z7 * x1 + z6 * x2 + z5 * x3;
  yqe(12) = z7 * y1 + z6 * y2 + z5 * y3;
  xqe(13) = ( x1 + x2 + x3 ) / 3.0;
  yqe(13) = ( y1 + y2 + y3 ) / 3.0;

  return
end
