function [ node_xy ] = xy_set ( nx, ny, node_num, xl, xr, yb, yt  )


  node_xy = zeros ( 2, node_num );

  for j = 1 : ny
    for i = 1 : nx 

      node_xy(1,i+(j-1)*(nx)) =    ...
        ( (  nx - i  ) * xl   ...
        + (          i - 1 ) * xr ) ...
        / (  nx     - 1 );

      node_xy(2,i+(j-1)*(nx)) =    ...
        ( (  ny - j ) * yb   ...
        + (          j - 1 ) * yt ) ...
        / ( ny     - 1);

    end
  end

  return
end
