function [ u, dudx, dudy ] = exact ( x, y )

%*****************************************************************************80
%
%% EXACT calculates the exact solution and its first derivatives.
%
%  Discussion:
%
%    The function specified here depends on the problem being
%    solved.  The user must make sure that EXACT and RHS are consistent.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    06 April 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, Y, the coordinates of a point
%    in the region, at which the exact solution is to be evaluated.
%
%    Output, real U, DUDX, DUDY, the value of
%    the exact solution U and its derivatives dUdX
%    and dUdY at the point (X,Y).
%
  u    =      sin ( pi * x ) * sin ( pi * y ) + x;
  dudx = pi * cos ( pi * x ) * sin ( pi * y ) + 1.0;
  dudy = pi * sin ( pi * x ) * cos ( pi * y );

  return
end