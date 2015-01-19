function [errL2,errH1] = fe_error2d(p, e, t, uh, u, u_x, u_y, k, t_point)

squareL2 = 0;
squareH1 = 0;
nt = size(t,2);

switch k
    case 1
        dof = 3;
        q = 1;
        for i = 1:nt
            [loc,glob,weights,DFinv] = quad_nodes(p(:,t(1,i)),p(:,t(2,i)),p(:,t(3,i)),q);
            nq = size(loc,2);
            for j = 1:nq
                uh_val = 0;
                grad_val = 0;
                for k = 1:dof
                    [v,grad] = P1ShapeFunction(k,loc(:,j));
                    uh_val = uh_val + uh(t(k,i))*v;
                    grad_val = grad_val + uh(t(k,i))*(DFinv*grad);
                end
                squareL2 = squareL2 + weights(j)*(u(glob(1,j),glob(2,j),t_point)-uh_val)^2;
                squareH1 = squareH1 + weights(j)*((u_x(glob(1,j),glob(2,j),t_point)-grad_val(1))^2 + (u_y(glob(1,j),glob(2,j),t_point)-grad_val(2))^2);
            end
        end
    case 2
        dof = 6;
        q = 4;
        np = size(p,2);
        edges = GetEdges(np,e,t);
        for i = 1:nt
            [loc,glob,weights,DFinv] = quad_nodes(p(:,t(1,i)),p(:,t(2,i)),p(:,t(3,i)),q);
            nq = size(loc,2);
            TtoE = MakeTtoE(t(:,i), edges);
            id = [t(:,i);np + TtoE(:)];
            for j = 1:nq
                uh_val = 0;
                grad_val = 0;
                for k = 1:dof
                    [v,grad] = P2ShapeFunction(k,loc(:,j));
                    uh_val = uh_val + uh(id(k))*v;
                    grad_val = grad_val + uh(id(k))*(DFinv*grad);
                end
                squareL2 = squareL2 + weights(j)*(u(glob(1,j),glob(2,j),t_point)-uh_val)^2;
                squareH1 = squareH1 + weights(j)*((u_x(glob(1,j),glob(2,j),t_point)-grad_val(1))^2 + (u_y(glob(1,j),glob(2,j),t_point)-grad_val(2))^2);
            end
        end
end

errL2 = sqrt(squareL2);
errH1 = sqrt(squareL2 + squareH1);