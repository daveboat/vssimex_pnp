function [dt_thresh,root_below,root_above,nogood,dt_close] = find_roots(j,epsilon,dt_end,root_tol,N,dt,dt_thresh,root_below,root_above,nogood,dt_close,v,dt_now,dt_old,phi_x_right,cp_ss,cm_ss,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)
    dt_now = .95*dt_end(j);
    dt_now = dt_close(j);

    
    % save on typing below (sheer laziness)
    cm = cm_ss;
    cp = cp_ss;
    test = 'quadratic';
    h = .01;
    
    [B_new,B_now,B_old] = auto_linearize_RHS(h,test,v,v,dt,dt,phi_x_right,phi_x_right,cp,cp,cp,cm,cm,cm,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
   
        
    root_old = 0;
    root_now = .5;
    ii=1;
    while abs(root_now) < 1 
        M_new = (3/2)*eye(2*N) - dt_now*B_new; % tridiagonal
        M_new = inv(M_new);
        M_now = 2*eye(2*N) + dt_now*B_now; % upper triangular
        M_old = (-1/2)*eye(2*N) + dt_now*B_old; % upper triangular
        A = [ M_new*M_now , M_new*M_old ; eye(2*N) , zeros(2*N,2*N)];
        [V,d,flag] = eigs(A,10);
        % did_they_converge(ii) = flag;
        % find and zero-out the eigenvalue that equals 1.
        a = abs(diag(d)); % sort the eigenvalues by magnitude
        I = find(abs(a-1) < 10^(-8)); % find the eigenvalue of magnitude 1
        d(I,I) = 0; % also set it to zero within the diagonal matrix
        % now find the largest-magnitude eigenvalue
        [a,I] = sort(abs(diag(d))); % sort the eigenvalues by magnitude    
        
        % update the roots and the times
        root_older = root_old;
        root_old = root_now;
        root_now = d(I(end),I(end));
        dt_older = dt_old;
        dt_old = dt_now;
        dt_now = dt_now + dt;
        
        ii = ii+1;
    end
    abs(root_now)-abs(root_old)
    ii
    % if we've exited this loop then we've found a value of dt which has its
    % largest root with magnitude greater than 1
    if abs(root_now)-abs(root_old) > root_tol
        nogood(j) = 1;
        display('fail')
        dt_close(j) = dt_older;
    else
        dt_thresh(j) = (dt_old+dt_now)/2;
        root_below(j) = root_old;
        root_above(j) = root_now;
    end
