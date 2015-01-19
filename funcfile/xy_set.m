function [ node_xy ] = xy_set ( nx, ny, node_num, xl, xr, yb, yt  )

%*****************************************************************************80
%
%% XY_SET sets the XY coordinates of the nodes.
%
%  Discussion:
%
%    The nodes are laid out in an evenly spaced grid, in the unit square.
%
%    The first node is at the origin.  More nodes are created to the
%    right until the value of X = 1 is reached, at which point
%    the next layer is generated starting back at X = 0, and an
%    increased value of Y.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 March 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer NX, NY, the number of elements in the X and
%    Y direction.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, real XL, XR, YB, YT, the X coordinates of
%    the left and right sides of the rectangle, and the Y coordinates
%    of the bottom and top of the rectangle.
%
%    Output, real NODE_XY(2,NODE_NUM), the nodes.
%
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
