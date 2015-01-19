function [ u, dudx, dudy ] = exact ( x, y )


  u    =      sin ( pi * x ) * sin ( pi * y ) + x;
  dudx = pi * cos ( pi * x ) * sin ( pi * y ) + 1.0;
  dudy = pi * sin ( pi * x ) * cos ( pi * y );

  return
end