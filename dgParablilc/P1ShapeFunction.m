function [v, grad] = P1ShapeFunction(k, xi)
  switch (k) 
   case 1 
    v = 1 - xi(1) - xi(2);
    grad = [-1; -1];
   case 2
    v = xi(1);
    grad = [1; 0];
   case 3
    v = xi(2);
    grad = [0; 1];
  end
  
