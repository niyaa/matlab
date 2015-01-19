function compare ( node_num, node_xy, f )


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