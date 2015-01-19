function [ wq, xq, yq ] = quad_a ( node_xy, element_node, ...
  element_num, node_num, element_order )


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
