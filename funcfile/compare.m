function compare ( node_num, node_xy, f )

%*****************************************************************************80
%
%% COMPARE compares the exact and computed solution at the nodes.
%
%  Discussion:
%
%    This is a rough comparison, done only at the nodes.  Such a pointwise
%    comparison is easy, because the value of the finite element
%    solution is exactly the value of the finite element coefficient
%    associated with that node.
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
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, real NODE_XY(2,NODE_NUM), the nodes.
%
%    Input, real F(NODE_NUM), the solution vector of the finite
%    element system.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'COMPARE:\n' );
  fprintf ( 1, '  Compare computed and exact solutions at the nodes.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '         X           Y          U           U\n' );
  fprintf ( 1, '                             computed     exact\n' );
  fprintf ( 1, '\n' );

  for node = 1 : node_num

    x = node_xy(1,node);
    y = node_xy(2,node);

    u = exact ( x, y );

    uh = f(node);

    fprintf ( 1, '%12f  %12f  %12f  %12f\n', x, y, uh, u );

  end

  return
end