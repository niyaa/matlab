function [ el2, eh1 ] = errors ( element_area, element_node, node_xy, ...
  f, element_num, element_order, nqe, node_num )

  el2 = 0.0E+00;
  eh1 = 0.0E+00;

  for element = 1 : element_num

    [ wqe, xqe, yqe ] = quad_e ( node_xy, element_node, element, ...
      element_num, element_order, node_num, nqe );

    for iq = 1 : nqe

      ar = element_area(element) * wqe(iq);
      x = xqe(iq);
      y = yqe(iq);

      uh = 0.0E+00;
      dudxh = 0.0E+00;
      dudyh = 0.0E+00;

      for in1 = 1 : element_order

        i = element_node(in1,element);

        [ bi, dbidx, dbidy ] = qbf ( x, y, element, in1, node_xy, ...
          element_node, element_num, element_order, node_num );

        uh    = uh    + bi    * f(i);
        dudxh = dudxh + dbidx * f(i);
        dudyh = dudyh + dbidy * f(i);

      end

      [ u, dudx, dudy ] = exact ( x, y );

      el2 = el2 + ( uh - u )^2 * ar;
      eh1 = eh1 + ( ( dudxh - dudx )^2 + ( dudyh - dudy )^2 ) * ar;

    end

  end

  el2 = sqrt ( el2 );
  eh1 = sqrt ( eh1 );

  fprintf ( 1, '\n' );
  fprintf ( 1, '*********************************************\n' );
  fprintf ( 1, '*                                           *\n' );
  fprintf ( 1, '*  ERRORS:                                  *\n' );
  fprintf ( 1, '*    L2 error =          %14f     *\n', el2 );
  fprintf ( 1, '*    H1-seminorm error = %14f     *\n', eh1 );
  fprintf ( 1, '*                                           *\n' );
  fprintf ( 1, '*********************************************\n' );

  return
end