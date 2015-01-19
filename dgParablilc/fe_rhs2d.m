function rhs = fe_rhs2d(p, e, t, f, bc, g, k, time_deg, time, u0, M)

% Calculates the right hand side for one time interval

np = size(p, 2);
nt = size(t, 2);

switch k % polynomial degree of spatial finite elements
    case 1
        dof_loc = 3;
        q = 2;
        shapefunction = @P1ShapeFunction;
        dof_glob = np;
    case 2
        dof_loc = 6;
        q = 4;
        shapefunction = @P2ShapeFunction;
        edges = GetEdges(np, e, t);
        TtoE = MakeTtoE(t, edges);
        dof_glob = np + size(edges, 2);
end

switch time_deg % polynomial degree of discretisation in time
    case 0
        rhs = zeros(dof_glob,1);
    case 1
        rhs = zeros(2*dof_glob,1); % RHS contains two blocks F1 and F2  
end

for i = 1:nt 
    p1 = p(:,(t(1,i))); 
    p2 = p(:,(t(2,i)));
    p3 = p(:,(t(3,i)));
    rhs_loc = rhs_element_matrix2d(p1, p2, p3, f, q, dof_loc, ...
                                 shapefunction, time_deg, time);  
    switch k
        case 1
            idx = t(1:3,i); 
        case 2
            idx = [t(1:3,i)' np+TtoE(1:3,i)']; 
    end                             
    switch time_deg
        case 0
            rhs(idx) = rhs(idx) + rhs_loc;
        case 1
            % block F1
            rhs(idx) = rhs(idx) + rhs_loc(1:dof_loc); 
            % block F2
            rhs(dof_glob+idx) = rhs(dof_glob+idx) + rhs_loc(1+dof_loc:2*dof_loc);
    end
end

switch time_deg
    case 0
        rhs = rhs + M*u0;
    case 1
        % block F1
        rhs(1:dof_glob) = rhs(1:dof_glob) + M*u0;
        % block F2
        rhs(1+dof_glob:2*dof_glob) = rhs(1+dof_glob:2*dof_glob)/diff(time);
end

%%% LOCAL FUNCTIONS %%%

function rhs_loc = rhs_element_matrix2d(p1, p2, p3, f, q, dof,...
                                        shapefunction, time_deg, time)

switch time_deg
    case 0
        rhs_loc = zeros(dof, 1);
        t_j = mean(time);
        delta_t = diff(time);
    case 1
        rhs_loc = zeros(2*dof, 1); % RHS contains two blocks
        t_j = mean(time);          % integration data for first block  
        delta_t = diff(time); 
        t_ref = [-1/sqrt(3) 1/sqrt(3)]; % integration data for second block
        t_points = 0.5 * (delta_t*t_ref + sum(time));
        w_t = 0.5 * delta_t;
end
                                        
[loc, glob, weights, DFinv] = quad_nodes(p1, p2, p3, q);

for j = 1:length(weights)
    x_j = glob(1, j);
    y_j = glob(2, j);
    w_j = weights(j);
    switch time_deg
        case 0
            f_j = f(x_j, y_j, t_j);
            for i = 1:dof
                [v_i, sgrad_i] = shapefunction(i, loc(:, j));
                rhs_loc(i) = rhs_loc(i) + w_j * v_i * delta_t * f_j;
            end
        case 1
            f_j = f(x_j, y_j, t_j);
            f_points = [f(x_j, y_j, t_points(1)), f(x_j, y_j, t_points(2))];
            for i = 1:dof
                [v_i, sgrad_i] = shapefunction(i, loc(:, j));
                % integration of first block
                rhs_loc(i) = rhs_loc(i) + w_j * v_i * delta_t * f_j; 
                % integration of second block
                rhs_loc(dof+i) = rhs_loc(dof+i) + w_j * v_i * w_t * ...
                                 sum((t_points-time(1)).*f_points);
            end
    end
end
 

