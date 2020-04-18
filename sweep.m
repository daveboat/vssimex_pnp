% Now try to compute that threshold value of dt...  At first I tried a
% naive approach and found it didn't work well enough.
%
% clear
% 
% load steady_state_data.mat
% root_tol = .1;
% 
% for j=1:length(EPSILON)
%     j
%     dt_now = .8*dt_end(j);
%     epsilon = EPSILON(j);
%     parameters
%     cm = cm_ss(:,j);
%     cp = cp_ss(:,j);
%     phi_x_right = phi_x_right_ss(j);
%     phi = phi_ss(:,j);
%     test = 'quadratic';
%     h = .01;
%     
%     [B_new,B_now,B_old] = auto_linearize_RHS(h,test,v(j),v(j),dt,dt,phi_x_right,phi_x_right,cp,cp,cp,cm,cm,cm,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
%     
%     dt = (2*dt_end(j)-.1*dt_end(j))/100;
%     
%     root_old = 0;
%     root_now = .5;
%     ii=1;
%     while abs(root_now) < 1
%         M_new = (3/2)*eye(2*N) - dt_now*B_new; % tridiagonal
%         M_new = inv(M_new);
%         M_now = 2*eye(2*N) + dt_now*B_now; % upper triangular
%         M_old = (-1/2)*eye(2*N) + dt_now*B_old; % upper triangular
%         A = [ M_new*M_now , M_new*M_old ; eye(2*N) , zeros(2*N,2*N)];
%         [V,d,flag] = eigs(A,10);
%         % did_they_converge(ii) = flag;
%         % find and zero-out the eigenvalue that equals 1.
%         a = abs(diag(d)); % sort the eigenvalues by magnitude
%         I = find(abs(a-1) < 10^(-8)); % find the eigenvalue of magnitude 1
%         d(I,I) = 0; % also set it to zero within the diagonal matrix
%         % now find the largest-magnitude eigenvalue
%         [a,I] = sort(abs(diag(d))); % sort the eigenvalues by magnitude
%         
%         % update the roots and the times
%         % root_older = root_old;
%         root_old = root_now;
%         root_now = d(I(end),I(end));
%         dt_old = dt_now;
%         dt_now = dt_now + dt;
%         
%         ii = ii+1;
%     end
%     dt_thresh(j) = (dt_old+dt_now)/2;
%     root_below(j) = root_old;
%     root_above(j) = root_now;
% end
% Let's see how well we approximated the critical time-step size.  Look
% at the following plots --- it's clear that the naive approach isn't doing
% too well.
% 
% figure(1)
% clf
% plot(EPSILON,abs(root_below),'r')
% hold on
% plot(EPSILON,abs(root_above),'r')
% title('deviation between root_{above} and root_{below}')
% xlabel('epsilon')
% % %
% figure(2)
% clf
% plot(EPSILON,(dt_thresh),'-')
% hold on
% plot(EPSILON,(dt_end),'r-')
% title('red: dt_{thresh} from adaptive timestepper, blue: from eigenvalues')
% xlabel('epsilon')
% ylabel('dt_{thresh}','FontSize',14)
% % %
% % 
% figure(3)
% plot(EPSILON,log10(abs(dt_thresh-dt_end)))

% pause

% What's going on?  When I try to compute dt_thresh from the linear
% stability analysis, I start with a dt_now that is below the dt_thresh from
% the adaptive time-stepper.  I compute the largest magnitude eigenvalue
% for that dt.  I then increase dt_now by dt and try again.  I keep doing this until
% I find a dt_old such that the largest magnitude eigenvalue (root_below) has
% magnitude less than 1 and dt_now = dt_old + dt such that the largest magnitude
% eigenvalue (root_above) has magnitude greater than 1.  I then define
% dt_thresh to be (dt_old+dt_now)/2 = dt_old + dt/2 .  The problem is: if I
% do all this with a fixed value of dt then I find that for relatively
% large values of epsilon, the difference between root_above and root_below
% is pretty small --- they both have magnitudes close to 1 and so averaging
% dt_old and dt_now isn't a bad idea.  But for smaller values of epsilon,
% there's a big gap between root_below and root_above (see figure 1) and
% this causes problems when looking at dt_thresh versus dt_end as a
% function of epsilon (see figures 2 and 3).
%
% So what I need to do to do things with smaller values of dt for small
% epsilon and larger values of dt for larger epsilon.  From Figure 1, I see
% that the gap is basically decreasing as epsilon decreases and so what
% I'll do is do the problem in reverse direction (start with the large
% values of epsilon and then go small) and I'll first do a run with a fixed
% dt and identify those values of epsilon for which dt was too large.  This
% is the "nogood" vector.  Now that I know which are the values of epsilon
% to focus on, I repeat this process with smaller and smaller values of dt
% until the nogood vector is an empty vector.  That's what's done below.

clear
load steady_state_data.mat
clear nogood
root_tol = .0001; % define how close i'd like abs(root_below) and abs(root_above) to be
% the smaller root_tol is, the longer it takes things to run.  For the
% figures in the article, I took root_tol = .00001.  I also had a finer array of 
% epsilon values than the range given here, especially close to the values
% where the corner and jump occur
%
% for j=1:length(Epsilon)
dt = .000001;
dt_thresh = zeros(size(EPSILON));
root_above = .5*ones(size(EPSILON));
root_below = zeros(size(EPSILON));
for j=length(EPSILON):-1:1    
    dt_now = .95*dt_end(j);
    epsilon = EPSILON(j);
    parameters
    cm = cm_ss(:,j);
    cp = cp_ss(:,j);
    phi_x_right = phi_x_right_ss(j);
    phi = phi_ss(:,j);
    test = 'quadratic';
    h = .01;
    
    [B_new,B_now,B_old] = auto_linearize_RHS(h,test,v(j),v(j),dt,dt,phi_x_right,phi_x_right,cp,cp,cp,cm,cm,cm,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    % dt = (2*dt_end(j)-.1*dt_end(j))/100;
    
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
    % if we've exited this loop then we've found a value of dt which has its
    % largest root with magnitude greater than 1
    if abs(root_now)-abs(root_old) > root_tol
        nogood(j) = 1;
        dt_close(j) = dt_older;
    else
        dt_thresh(j) = (dt_old+dt_now)/2;
        root_below(j) = root_old;
        root_above(j) = root_now;
    end
    j
end
% now that nogood is defined, I start decreasing dt.  I can look at these
% figures to see how things look...
% figure(1)
% clf
% plot(EPSILON((length(nogood)+1):end),dt_thresh((length(nogood)+1):end))
% hold on
% plot(EPSILON((length(nogood)+1):end),dt_end((length(nogood)+1):end),'ro')
% figure(2)
% plot(EPSILON((length(nogood)+1):end),log10(abs(root_below((length(nogood)+1):end)-root_above((length(nogood)+1):end))))

while length(nogood) > 0
    dt = dt/10;
    nogood = zeros(size(nogood));
    length(nogood)
    for j=length(nogood):-1:1
        epsilon = EPSILON(j);
        j
        parameters
        [dt_thresh,root_below,root_above,nogood,dt_close] = find_roots(j,epsilon,dt_end,root_tol,N,dt,dt_thresh,root_below,root_above,nogood,dt_close,v(j),dt_now,dt_old,phi_x_right_ss(j),cp_ss(:,j),cm_ss(:,j),bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    end
    nogood = find(nogood==1);
end

save steady_state_data.mat

% % Now lets repeat those figures and see how they look!
% figure(1)
% clf
% plot(abs(root_below))
% hold on
% plot(abs(root_above))
% title('root_above, root_below uniform-ish deviation for small epsilons')
% xlabel('epsilon')
% 
% 
% figure(2)
% clf
% plot(EPSILON,(dt_thresh),'.-')
% hold on
% plot(EPSILON,(dt_end),'r.-')
% title('red: dt_{thresh} from adaptive timestepper, blue: from eigenvalues')
% xlabel('epsilon')
% ylabel('dt_{thresh}','FontSize',14)

% % look at it on log10 scale for a better look...
% figure(3)
% clf
% plot(Epsilon,log10(dt_thresh),'o-')
% hold on
% plot(Epsilon,log10(dt_end),'ro-')
% title('red: dt_{thresh} from adaptive timestepper, blue: from eigenvalues')
% xlabel('epsilon')
% ylabel('log10(dt_{thresh})','FontSize',14)
