
function element_node = grid_t6 ( nx, ny, element_order, element_num )


  element = 0;

  for j = 1: ny - 1
    for i = 1 : nx - 1

      sw = ( j - 1 ) * ( nx  ) + i;
       se=sw+1;

      nw  = sw +  nx ;
      ne=nw+1;
    

      element = element + 1;
      element_node(1,element) = sw;
      element_node(2,element) = se;
      element_node(3,element) = nw;


      element = element + 1;
      element_node(1,element) = ne;
      element_node(2,element) = nw;
      element_node(3,element) = se;
 

    end
  end
  return
end