function value = rhs ( x, y )

%*****************************************************************************80
%
%% RHS gives the right-hand side of the differential equation.
%
%  Discussion:
%
%    The function specified here depends on the problem being
%    solved.  This is one of the routines that a user will
%    normally want to change.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 February 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, Y, the coordinates of a point
%    in the region, at which the right hand side of the
%    differential equation is to be evaluated.
%
%    Output, real VALUE, the value of the right
%    hand side of the differential equation at (X,Y).
%
  value = 2.0 * pi * pi * sin ( pi * x ) * sin ( pi * y );

  return
end