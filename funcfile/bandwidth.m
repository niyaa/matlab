function nhba = bandwidth ( element_order, element_num, element_node )

%*****************************************************************************80
%
%% BANDWIDTH determines the bandwidth of the coefficient matrix.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 April 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer ELEMENT_ORDER, the number of local nodes per element.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
%    ELEMENT_NODE(I,J) is the global index of local node I in element J.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, integer NHBA, the half bandwidth of the matrix.
%
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