function u = full_time_solution(p, e, t, alpha, gamma,...
                                f, bc, g, k, time_deg, time_vec, u0)

[M_gamma, M, S] = fe_matrices2d(p, e, t, alpha, gamma, k);
A = M_gamma + S;
ne = size(e,2);
np = size(p,2);
dof = length(u0);

for time_id = 1:length(time_vec)-1
    time = time_vec(time_id:time_id+1);
    delta_t = diff(time);
    rhs = fe_rhs2d(p, e, t, f, bc, g, k, time_deg, time, u0, M);
    switch time_deg
        case 0
            % We have to solve a system
            %
            %  P * u = F
            %
            % which yields the solution u
            P = M + delta_t*A;
            [P,rhs] = apply_bc(P, rhs, p, e, t, bc, g, k, time_deg, time);
            u = P\rhs;
        case 1
            % We have to solve a system
            %
            %    P1 Q1  *  u0  =  F1
            %    P2 Q2     u1     F2
            %
            % which yields the solution u = u0 + u1
            P1 = M + delta_t*A;
            Q1 = M + 0.5*delta_t*A;
            P2 = 0.5*delta_t*A;
            Q2 = 0.5*M + (delta_t/3)*A;
            P = [P1 Q1; P2 Q2];
            [P,rhs] = apply_bc(P, rhs, p, e, t, bc, g, k, time_deg, time, u0);
            u_pre = P\rhs;
            u = u_pre(1:dof) + u_pre(1+dof:2*dof);
    end
    u0 = u;
end