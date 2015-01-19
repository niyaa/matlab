function [ wq, xq, yq ] = quad_a ( node_xy, element_node, ...
  element_num, node_num, element_order )

%*****************************************************************************80
%
%% QUAD_A sets the quadrature rule for assembly.
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
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
%
%    Output, real WQ(3), quadrature weights.
%
%    Output, real XQ(3,ELEMENT_NUM), YQ(3,ELEMENT_NUM), the X and Y
%    coordinates of the quadrature points in each element.
%
  wq = zeros ( 3, 1 );
  xq = zeros ( 3, element_num );
  yq = zeros ( 3, element_num );

  wq(1) = 1.0 / 3.0;
  wq(2) = wq(1);
  wq(3) = wq(1);

  for element = 1 : element_num

    ip1 = element_node(1,element);
    ip2 = element_node(2,element);
    ip3 = element_node(3,element);

    x1 = node_xy(1,ip1);
    x2 = node_xy(1,ip2);
    x3 = node_xy(1,ip3);

    y1 = node_xy(2,ip1);
    y2 = node_xy(2,ip2);
    y3 = node_xy(2,ip3);

    xq(1,element) = 0.5E+00 * ( x1 + x2 );
    xq(2,element) = 0.5E+00 * ( x2 + x3 );
    xq(3,element) = 0.5E+00 * ( x1 + x3 );

    yq(1,element) = 0.5E+00 * ( y1 + y2 );
    yq(2,element) = 0.5E+00 * ( y2 + y3 );
    yq(3,element) = 0.5E+00 * ( y1 + y3 );

  end

  return
end
