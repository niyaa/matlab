function nhba = bandwidth ( element_order, element_num, element_node )


  nhba = 0;
  for element = 1 : element_num
    for iln = 1 : element_order
      i = element_node(iln,element);
      if ( 0 < i )
        for jln = 1 : element_order
          j = element_node(jln,element);
          nhba = max ( nhba, j - i );
        end
      end
    end
  end

  return
end