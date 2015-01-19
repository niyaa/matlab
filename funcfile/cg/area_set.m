function element_area = area_set ( node_xy, element_num, element_node )


  element_area = zeros ( element_num, 1 );

  for element = 1 : element_num

    i1 = element_node(1,element);
    x1 = node_xy(1,i1);
    y1 = node_xy(2,i1);

    i2 = element_node(2,element);
    x2 = node_xy(1,i2);
    y2 = node_xy(2,i2);

    i3 = element_node(3,element);
    x3 = node_xy(1,i3);
    y3 = node_xy(2,i3);

    element_area(element) = 0.5 * abs ...
      ( y1 * ( x2 - x3 ) ...
      + y2 * ( x3 - x1 ) ...
      + y3 * ( x1 - x2 ) );

  end

  return
end