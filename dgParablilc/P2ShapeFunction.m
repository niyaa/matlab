function [v, grad] = P2ShapeFunction(k, xi)
  switch (k) 
   case  1
    v = (1 - xi(1) - xi(2))*(1 - 2*xi(1) - 2*xi(2));
    grad = [-3 + 4*xi(1) + 4*xi(2); -3 + 4*xi(2) + 4*xi(1)];
   case  2
    v = xi(1)*(2*xi(1) - 1);
    grad = [4*xi(1) - 1; 0];
   case  3
    v = xi(2)*(2*xi(2) - 1);	
    grad = [0; 4*xi(2) - 1];
   case  4
    v = 4*xi(1)*(1 - xi(1) - xi(2));
    grad = [4 - 8*xi(1) - 4*xi(2); -4*xi(1)];
   case  5
    v = 4*xi(1)*xi(2);
    grad = [4*xi(2); 4*xi(1)];
   case  6
    v = 4*xi(2)*(1 - xi(1) - xi(2));
    grad = [-4*xi(2); 4 - 8*xi(2) - 4*xi(1)];
  end
  
