
function element_node = grid_t6 ( nx, ny, nnodes, element_num )


  element = 0;

  for j = 1: ny - 1
    for i = 1 : nx - 1

      sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 1;
      w  = sw + 1;
      nw = sw + 2;

      s  = sw + 2 * nx - 1;
      c  = s  + 1;
      n  = s  + 2;

      se = s  + 2 * nx - 1;
      e  = se + 1;
      ne = se + 2;

      element = element + 1;
      element_node(1,element) = sw;
      element_node(2,element) = se;
      element_node(3,element) = nw;
      element_node(4,element) = s;
      element_node(5,element) = c;
      element_node(6,element) = w;

      element = element + 1;
      element_node(1,element) = ne;
      element_node(2,element) = nw;
      element_node(3,element) = se;
      element_node(4,element) = n;
      element_node(5,element) = c;
      element_node(6,element) = e;

    end
  end
  return
end