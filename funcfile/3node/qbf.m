
function [ b, dbdx, dbdy ] = qbf ( x, y, element, inode, node_xy, ...
  element_node, element_num, element_order, node_num )


%
  xn(1:3) = node_xy(1,element_node(1:3,element));
  yn(1:3) = node_xy(2,element_node(1:3,element));
%
%  Determine the (R,S) coordinates corresponding to (X,Y).
%
%  What is happening here is that we are solving the linear system:
%
%    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
%    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
%
%  by computing the inverse of the coefficient matrix and multiplying
%  it by the right hand side to get R and S.
%
%  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
%  for R and S.
%
  det =   ( xn(2) - xn(1) ) * ( yn(3) - yn(1) ) ...
        - ( xn(3) - xn(1) ) * ( yn(2) - yn(1) );

  r = ( ( yn(3) - yn(1) ) * ( x     - xn(1) ) ...
      + ( xn(1) - xn(3) ) * ( y     - yn(1) ) ) / det;

  drdx = ( yn(3) - yn(1) ) / det;
  drdy = ( xn(1) - xn(3) ) / det;

  s = ( ( yn(1) - yn(2) ) * ( x     - xn(1) ) ...
      + ( xn(2) - xn(1) ) * ( y     - yn(1) ) ) / det;

  dsdx = ( yn(1) - yn(2) ) / det;
  dsdy = ( xn(2) - xn(1) ) / det;
%
%  The basis functions can now be evaluated in terms of the
%  reference coordinates R and S.  It's also easy to determine
%  the values of the derivatives with respect to R and S.
%
  if ( inode == 1 )

    b    =   1;
    dbdr = 0.0;
    dbds = 0.0;

  elseif ( inode == 2 )

    b    =   r;
    dbdr = 1.0 ;
    dbds =   0.0;

  elseif ( inode == 3 )

    b    =   s ;
    dbdr =   0.0;
    dbds = 1.0 ;

  

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'QBF - Fatal error!\n' );
    fprintf ( 1, '  Request for local basis function INODE = %d\n', inode );
    error ( 'QBF - Fatal error!' );

  end
%
%  We need to convert the derivative information from (R(X,Y),S(X,Y))
%  to (X,Y) using the chain rule.
%
  dbdx = dbdr * drdx + dbds * dsdx;
  dbdy = dbdr * drdy + dbds * dsdy;

  return
end


