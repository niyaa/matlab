function [ node_xy ] = xy_set ( nx, ny, node_num, xl, xr, yb, yt  )


  node_xy = zeros ( 2, node_num );

  for j = 1 : 2*ny-1
    for i = 1 : 2*nx - 1

      node_xy(1,i+(j-1)*(2*nx-1)) =    ...
        ( ( 2 * nx - i - 1 ) * xl   ...
        + (          i - 1 ) * xr ) ...
        / ( 2 * nx     - 2 );

      node_xy(2,i+(j-1)*(2*nx-1)) =    ...
        ( ( 2 * ny - j - 1 ) * yb   ...
        + (          j - 1 ) * yt ) ...
        / ( 2 * ny     - 2 );

    end
  end

  return
end
