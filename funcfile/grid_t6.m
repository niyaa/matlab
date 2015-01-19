
function element_node = grid_t6 ( nx, ny, nnodes, element_num )

%*****************************************************************************80
%
%% GRID_T6 produces a grid of pairs of 6 node triangles.
%
%  Example:
%
%    Input:
%
%      NX = 4, NY = 3
%
%    Output:
%
%      ELEMENT_NODE =
%         1,  3, 15,  2,  9,  8;
%        17, 15,  3, 16,  9, 10;
%         3,  5, 17,  4, 11, 10;
%        19, 17,  5, 18, 11, 12;
%         5,  7, 19,  6, 13, 12;
%        21, 19,  7, 20, 13, 14;
%        15, 17, 29, 16, 23, 22;
%        31, 29, 17, 30, 23, 24;
%        17, 19, 31, 18, 25, 24;
%        33, 31, 19, 32, 25, 26;
%        19, 21, 33, 20, 27, 26;
%        35, 33, 21, 34, 27, 28.
%
%  Diagram:
%
%   29-30-31-32-33-34-35
%    |\ 8  |\10  |\12  |
%    | \   | \   | \   |
%   22 23 24 25 26 27 28
%    |   \ |   \ |   \ |
%    |  7 \|  9 \| 11 \|
%   15-16-17-18-19-20-21
%    |\ 2  |\ 4  |\ 6  |
%    | \   | \   | \   |
%    8  9 10 11 12 13 14
%    |   \ |   \ |   \ |
%    |  1 \|  3 \|  5 \|
%    1--2--3--4--5--6--7
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
%    Input, integer NX, NY, controls the number of elements along the
%    X and Y directions.  The number of elements will be
%    2 * ( NX - 1 ) * ( NY - 1 ).
%
%    Input, integer NNODES, the number of local nodes per element.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Output, integer ELEMENT_NODE(NNODES,ELEMENT_NUM);
%    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
%
  element = 0;

  for j = 1: ny - 1
    for i = 1 : nx - 1

      sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 1;
      w  = sw + 1;
      nw = sw + 2;

      s  = sw + 2 * nx - 1;
      c  = s  + 1;
      n  = s  + 2;

      se = s  + 2 * nx - 1;
      e  = se + 1;
      ne = se + 2;

      element = element + 1;
      element_node(1,element) = sw;
      element_node(2,element) = se;
      element_node(3,element) = nw;
      element_node(4,element) = s;
      element_node(5,element) = c;
      element_node(6,element) = w;

      element = element + 1;
      element_node(1,element) = ne;
      element_node(2,element) = nw;
      element_node(3,element) = se;
      element_node(4,element) = n;
      element_node(5,element) = c;
      element_node(6,element) = e;

    end
  end
  return
end