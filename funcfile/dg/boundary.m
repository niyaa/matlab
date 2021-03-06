function [ a, f ] = boundary ( nx, ny, node_num, node_xy, ib, a, f )


  node = 0;

  for row = 1 : 2 * ny - 1

    for col = 1 : 2 * nx - 1

      node = node + 1;

      if ( row == 1 || row == 2 * ny - 1 || col == 1 || col == 2 * nx - 1 )

        x = node_xy(1,node);
        y = node_xy(2,node);
        u = exact ( x, y );

        jlo = max ( node - ib, 1 );
        jhi = min ( node + ib, node_num );

        for j = jlo : jhi
          a(node,j) = 0.0;
        end

        a(node,node) = 1.0;

        f(node) = u;

      end

    end
  end

  return
end