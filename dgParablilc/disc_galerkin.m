function disc_galerkin()

% displays the error decay for the numeric solution of the 2D Equation
%
% du/dt - div(alpha*grad(u)) + gamma*u = f
%
% by the Discontinuous Galerkin Method for either fixed h or fixed delta_t
% in the domain (x,y,t) \in [-1,1] x [-1,1] x [0,1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fixed = 'h';                        % fix either h or t
k = 1;                              % polynomial degree of finite elements in space: 1,2
time_deg = 1;                       % polynomial degree of finite elements in time: 0,1
alpha = @(x,y)[1 0; 0 1];           % PDE parameters
gamma = @(x,y)0;               
f = @(x,y,t)3*t^2*(x^2+y^2)-4*t^3;        
u_exact = @(x,y,t)(x^2+y^2)*t^3;    % exact solution
u_x = @(x,y,t)2*x*t^3;              % partial derivate wrt. x
u_y = @(x,y,t)2*y*t^3;              % partial derivate wrt. y
bc ='dddd';                         % boundary conditions (d: Dirichlet, n: Neumann)
g = @get_g;                         % boundary values (see below)
error = 'L2';                       % error to be plotted: L2,H1,EK (EK = L-infinity norm)

if fixed == 'h'
    h = 0.025;                      % fix h
    max_steps = 5;                  % maximal number of refinement steps for t
    ot_order = 2;                   % expected approximation order
elseif fixed == 't'
    delta_t = 0.01;                 % fix delta_t
    max_steps = 4;                  % maximal number of refinement steps for h
    oh_order = 1;                   % expected approximation order
end

errL2_vec = [];
errH1_vec = [];
errEK_vec = [];
disp('Discontinuous Galerkin Method');

% Decrease delta_t for a fixed h
if fixed == 'h'
    
    [p,e,t] = triangulate([-1 1 1 -1; -1 -1 1 1], h, 3);                    
    t_max = 1;
    np = size(p,2);
    if k == 1
        dof_glob = np;
    elseif k == 2
        edges = GetEdges(np, e, t);
        dof_glob = np + size(edges, 2);
    end    
    u0 = zeros(dof_glob,1);
    u_final = zeros(dof_glob,1);
    for i = 1:np
        u0(i) = u_exact(p(1,i),p(2,i),0);
        u_final(i) = u_exact(p(1,i),p(2,i),t_max);
    end
    if k == 2
        p_mean = 0.5*(p(:,edges(1,:))+p(:,edges(2,:)));
        for i = 1:size(edges, 2)
            u0(np+i) = u_exact(p_mean(1,i),p_mean(2,i),0);
            u_final(np+i) = u_exact(p_mean(1,i),p_mean(2,i),t_max);
        end
    end
    
    delta_t_vec = [];
    o_t_vec = [];

    for t_exp = 0:max_steps
        delta_t = 1/2^t_exp;
        disp(['delta_t = ',num2str(delta_t)]);
        time_vec = 0:delta_t:t_max;
        u = full_time_solution(p, e, t, alpha, gamma,...
                                f, bc, g, k, time_deg, time_vec, u0);
        [errL2,errH1] = fe_error2d(p, e, t, u, u_exact, u_x, u_y, k, t_max);                                               
        errEK = norm(u_final - u);
        errL2_vec = [errL2_vec, errL2];
        errH1_vec = [errH1_vec, errH1];
        errEK_vec = [errEK_vec, errEK];
        delta_t_vec = [delta_t_vec, delta_t];
        o_t_vec = [o_t_vec, 0.1 * delta_t^ot_order];
    end
    
    
% Decrease h for a fixed delta_t    
elseif fixed == 't'
    
    t_max = 1;
    time_vec = 0:delta_t:t_max;
 
    h_vec = [];
    o_h_vec = [];
    
    for h_exp = 0:max_steps
        h = 1/2^h_exp;
        disp(['h = ',num2str(h)]);
        [p,e,t] = triangulate([-1 1 1 -1; -1 -1 1 1], h, 3);
        
        np = size(p,2);
        if k == 1
            dof_glob = np;
        elseif k == 2
            edges = GetEdges(np, e, t);
            dof_glob = np + size(edges, 2);
        end    
        u0 = zeros(dof_glob,1);
        u_final = zeros(dof_glob,1);
        for i = 1:np
            u0(i) = u_exact(p(1,i),p(2,i),0);
            u_final(i) = u_exact(p(1,i),p(2,i),t_max);
        end
        if k == 2
            p_mean = 0.5*(p(:,edges(1,:))+p(:,edges(2,:)));
            for i = 1:size(edges, 2)
                u0(np+i) = u_exact(p_mean(1,i),p_mean(2,i),0);
                u_final(np+i) = u_exact(p_mean(1,i),p_mean(2,i),t_max);
            end
        end
        
        u = full_time_solution(p, e, t, alpha, gamma,...
                                    f, bc, g, k, time_deg, time_vec, u0);
        [errL2,errH1] = fe_error2d(p, e, t, u, u_exact, u_x, u_y, k, t_max);                                               
        errEK = norm(u_final - u);
        errL2_vec = [errL2_vec, errL2];
        errH1_vec = [errH1_vec, errH1];
        errEK_vec = [errEK_vec, errEK];
        h_vec = [h_vec, h];
        o_h_vec = [o_h_vec, 0.1 * h^oh_order];
    end
     
end


%%% Plot Options %%%

figure;

if fixed == 'h'
    
    loglog(delta_t_vec, o_t_vec, 'g--','LineWidth',3);
    set(gca,'XDir','reverse');
    set(gca,'FontSize',18)
    xlabel('\Deltat');
    title(['h = ', num2str(h),]); 
    hold on;
    if error == 'L2'
        loglog(delta_t_vec, errL2_vec, '-ro','LineWidth',3,'MarkerFaceColor','r','MarkerSize',5);
        ylabel('L^2 error');
    elseif error == 'H1'
        loglog(delta_t_vec, errH1_vec, '-ro','LineWidth',3,'MarkerFaceColor','r','MarkerSize',5);
        ylabel('H^1 error');
    elseif error == 'EK'
        loglog(delta_t_vec, errEK_vec, '-ro','LineWidth',3,'MarkerFaceColor','r','MarkerSize',5);
        ylabel('Eucl. error');
    end
    legend(['O(\Deltat^',num2str(ot_order),')'],'Disc-Galerkin','Location','SouthWest');
    hold off;
    
elseif fixed == 't'
    
    loglog(h_vec, o_h_vec, 'g--','LineWidth',3);
    set(gca,'XDir','reverse');
    set(gca,'FontSize',18)
    xlabel('h');
    title(['\Deltat = ', num2str(delta_t),]); 
    hold on;
    if error == 'L2'
        loglog(h_vec, errL2_vec, '-ro','LineWidth',3,'MarkerFaceColor','r','MarkerSize',5);
        ylabel('L^2 error');
    elseif error == 'H1'
        loglog(h_vec, errH1_vec, '-ro','LineWidth',3,'MarkerFaceColor','r','MarkerSize',5);
        ylabel('H^1 error');
    elseif error == 'EK'
        loglog(h_vec, errEK_vec, '-ro','LineWidth',3,'MarkerFaceColor','r','MarkerSize',5);
        ylabel('Eucl. error');
    end
    legend( ['O(h^',num2str(oh_order),')'],'Disc-Galerkin','Location','SouthWest');
    hold off;
        
end


%%% LOCAL FUNCTIONS %%%

function val = get_g(x, y, t, bs)

switch bs
    case {1,2,3,4}
        val = (x^2+y^2)*t^3;; % Dirichlet boundary conditions
end        