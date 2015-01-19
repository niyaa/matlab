function rhs_loc = fe_surfaceint(p1, p2, g, bs, k, time_deg, time)

% Calculates Neumann boundary conditions

switch time_deg
    case 0
        delta_t = diff(time);
        w_t = delta_t;
        switch k
            case 1
                phi = 0.5;
                p_x = 0.5*(p1+p2);
                w = norm(p1-p2);
                rhs_loc = ones(1,2) * w_t * w * phi * g(p_x(1),p_x(2),time(2),bs);
            case 2
                x = 0.5*[1-1/sqrt(3), 1+1/sqrt(3)];
                phi1 = 2*x.^2-3*x.+1;
                phi2 = 2*x.^2-x.;
                phi3 = -4*x.^2+4*x.;
                p_x1 = p1+x(1)*(p2-p1);
                p_x2 = p1+x(2)*(p2-p1);
                w = 0.5*norm(p1-p2);
                rhs_loc(1) = w_t * w * (phi1(1)*g(p_x1(1),p_x1(2),time(2),bs) + phi1(2)*g(p_x2(1),p_x2(2),time(2),bs));
                rhs_loc(2) = w_t * w * (phi2(1)*g(p_x1(1),p_x1(2),time(2),bs) + phi2(2)*g(p_x2(1),p_x2(2),time(2),bs));
                rhs_loc(3) = w_t * w * (phi3(1)*g(p_x1(1),p_x1(2),time(2),bs) + phi3(2)*g(p_x2(1),p_x2(2),time(2),bs));
        end
        
    case 1
        delta_t = diff(time); 
        t_ref = [-1/sqrt(3) 1/sqrt(3)];
        t_points = 0.5 * (delta_t*t_ref + sum(time));
        w_t = 0.5 * delta_t;
        switch k
            case 1
                phi = 0.5;
                p_x = 0.5*(p1+p2);
                w = norm(p1-p2);
                rhs_loc = ones(1,2) * w_t * w * phi * (g(p_x(1),p_x(2),t_points(1),bs) + g(p_x(1),p_x(2),t_points(2),bs);
            case 2
                x = 0.5*[1-1/sqrt(3), 1+1/sqrt(3)];
                phi1 = 2*x.^2-3*x.+1;
                phi2 = 2*x.^2-x.;
                phi3 = -4*x.^2+4*x.;
                p_x1 = p1+x(1)*(p2-p1);
                p_x2 = p1+x(2)*(p2-p1);
                w = 0.5*norm(p1-p2);
                rhs_loc(1) = w_t * w * ((phi1(1)*g(p_x1(1),p_x1(2),t_points(1),bs) + phi1(2)*g(p_x2(1),p_x2(2),t_points(1),bs)) + ...
                                        (phi1(1)*g(p_x1(1),p_x1(2),t_points(2),bs) + phi1(2)*g(p_x2(1),p_x2(2),t_points(2),bs));
                rhs_loc(2) = w_t * w * ((phi2(1)*g(p_x1(1),p_x1(2),t_points(1),bs) + phi2(2)*g(p_x2(1),p_x2(2),t_points(1),bs)) + ...
                                        (phi2(1)*g(p_x1(1),p_x1(2),t_points(2),bs) + phi2(2)*g(p_x2(1),p_x2(2),t_points(2),bs));
                rhs_loc(3) = w_t * w * ((phi3(1)*g(p_x1(1),p_x1(2),t_points(1),bs) + phi3(2)*g(p_x2(1),p_x2(2),t_points(1),bs)) + ...
                                        (phi3(1)*g(p_x1(1),p_x1(2),t_points(2),bs) + phi3(2)*g(p_x2(1),p_x2(2),t_points(2),bs));
        end
end