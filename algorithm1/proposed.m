function [theta_es,z_es,err]=proposed(dict,Y,X,Nx,Nt,Nr,Ny,Rth,Y_vec,Gaoa,Gaod,grid_aoa,grid_aod,L)
max_outer_iter=3;
max_inner_iter=500;
Snum1 = L;
Snum2 = L;
theta = zeros(2,0);
z_old = zeros(1,0);
for outer_iter = 1:max_outer_iter
    At=exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
    Ar=exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
    R = Y-Ar*diag(z_old)*At'*X;
    r=vec(R);
    Rnorm = norm(R,'fro');
    if ((outer_iter>1) && (Rnorm < Rth))
        break;
        % check if the residue is small enough
        % if so, break and return to the final result
        % else, maybe some paths are missed, such missing paths cannot be found by gradient descend
        %       then do SVD on the residue and find the missing paths and
        %       do IR-based channel estimation again
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U,~,V] = svd(R);
    % SVD based preconditioning
    % find the missed theta's in the previous iteration
    % add to the theta-list of the current iteration
    theta_update = zeros(2,Snum1);
%     for i = 1:Snum1
%         u = U(:,i);
%         [~,ui] = max(fft(u));
%         theta_update(2, i) = (Nr-ui+1)/Nr;
%         v = V(:,i);
%         [~,vi] = max(fft(X*v));
%         theta_update(1, i) = (Nt-vi+1)/Nt;
%     end
    
    for i=1:Snum1
        product = dict'*r;%传感矩阵A各列与残差的内积
        [val,pos] = max(abs(product));%找到最大内积绝对值，即与残差最相关的列
        Pos_t= floor((pos-1)/Gaoa)+1;%存储这一列的序号 先t再r
        theta_update(1, i) =grid_aod(Pos_t);
        Pos_r = mod(pos-1,Gaoa)+1;%存储这一列的序号
        theta_update(2, i) =grid_aoa(Pos_r);
    end

    theta = [theta theta_update];
    epsilon = 1;
    z_new=[z_old;ones(Snum1,1)];
    z_old=[z_old;zeros(Snum1,1)];

    index_amp=1:numel(z_old);

    stepsize_old = 1;
    for inner_itr=1:max_inner_iter
        if epsilon>1e-8 && norm(z_old-z_new)<epsilon^0.5
            epsilon=epsilon/sqrt(10);
        end
        z_old=z_new;
 
        dd=1./(abs(z_old).^2+epsilon);
        D=diag(dd);
   
        lambda = 10;
        % pruning and lambda update
        if epsilon<1e-3
            L_index_amp0 = length(index_amp);
            index_amp = 1:L_index_amp0;
            threshold=0.005;
            if (numel(z_new) > Snum2)
                z_sort = sort(abs(z_new),'descend');
                threshold = max(z_sort(Snum2+1),0.005);
            end
            index_t=find(abs(z_new)>threshold);
            if ~isempty(index_t)
                index_amp=index_t;
            end
            % prune the small components of z
            D=D(index_amp,index_amp);
            z_old = z_old(index_amp);
            theta = theta(:,index_amp);
            L_index_amp2 = numel(index_amp);
            At=exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
            Ar=exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
            R = Y-Ar*diag(z_old)*At'*X;
            Rnorm = norm(R,'fro');
            lambda = max( 1*(Rnorm^2),1e-8);
            % update the weight parameter lambda
        end

        theta_new = theta;
        L_new = length(index_amp);
        %% gradient descend
        dtheta = zeros(2, L_new);
        At = exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
        Ar = exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
        At_multiply_X = At'*X;
        At_multiply_X2 = (At_multiply_X*At_multiply_X').';
        Ar_multiply_W = Ar';
        Ar_multiply_W2 = Ar_multiply_W*Ar_multiply_W';
        sigma_ky = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky = sigma_ky + (Ar_multiply_W*Y(:,p)).*conj(At_multiply_X(:,p));
        end
        sigma_kk = At_multiply_X2.*Ar_multiply_W2;
        inv_dkk_multiply_sigma_ky = (D/lambda + sigma_kk)\(sigma_ky);
        pAt = diag(-2i*pi*(0:1:Nt-1))*exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
        pAr = diag(-2i*pi*(0:1:Nr-1))*exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
        pAt_multiply_X = pAt'*X;
        pAr_multiply_W = pAr';
        for i = 1:L_new
            p1 = zeros(1,L_new);
            p2 = zeros(1,L_new);
            for ii = 1:Nx
                p1(1,i) = p1(1,i) + Y(:,ii)'*Ar_multiply_W(i,:)'*pAt_multiply_X(i,ii);
                p2(1,i) = p2(1,i) + Y(:,ii)'*pAr_multiply_W(i,:)'*At_multiply_X(i,ii);
            end
            p30 = zeros(L_new,L_new);
            p40 = zeros(L_new,L_new);
            p30(:,i) = At_multiply_X*X'*pAt(:,i);
            p30(i,:) = p30(i,:) + p30(:,i)';
            p40(:,i) = Ar_multiply_W*pAr(:,i);
            p40(i,:) = p40(i,:) + p40(:,i)';
            p3 = p30.'.*Ar_multiply_W2;
            p4 = At_multiply_X2.*p40;

            df_dtheta_r_l = -p2*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p2'+inv_dkk_multiply_sigma_ky'*p4*inv_dkk_multiply_sigma_ky;
            df_dtheta_t_l = -p1*inv_dkk_multiply_sigma_ky-inv_dkk_multiply_sigma_ky'*p1'+inv_dkk_multiply_sigma_ky'*p3*inv_dkk_multiply_sigma_ky;

            dtheta(1, i) = real(df_dtheta_t_l);
            dtheta(2, i) = real(df_dtheta_r_l);
        end
        stepsize = stepsize_old*4;
        theta_new = theta - stepsize*dtheta;
        theta_new = mod(theta_new, 1);

        At_new = exp(-2i*pi*(0:1:Nt-1)'*theta_new(1,:));
        Ar_new = exp(-2i*pi*(0:1:Nr-1)'*theta_new(2,:));
        At_multiply_X_new = At_new'*X;
        At_multiply_X2_new = (At_multiply_X_new*At_multiply_X_new').';
        Ar_multiply_W_new = Ar_new';
        Ar_multiply_W2_new = Ar_multiply_W_new*Ar_multiply_W_new';
        sigma_ky_new = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky_new = sigma_ky_new + (Ar_multiply_W_new*Y(:,p)).*conj(At_multiply_X_new(:,p));
        end
        sigma_kk_new = At_multiply_X2_new.*Ar_multiply_W2_new;

        func_val = -sigma_ky_new'*((D/lambda + sigma_kk_new)\sigma_ky_new);
        sur_val = -sigma_ky'*inv_dkk_multiply_sigma_ky;

        maxit=0;
        while func_val>sur_val && maxit<50
            stepsize = 0.1*stepsize;
            theta_new = theta - stepsize*(dtheta/norm(dtheta,'fro'));
            theta_new = mod(theta_new, 1);

            At_new = exp(-2i*pi*(0:1:Nt-1)'*theta_new(1,:));
            Ar_new = exp(-2i*pi*(0:1:Nr-1)'*theta_new(2,:));
            At_multiply_X_new = At_new'*X;
            At_multiply_X2_new = (At_multiply_X_new*At_multiply_X_new').';
            Ar_multiply_W_new = Ar_new';
            Ar_multiply_W2_new = Ar_multiply_W_new*Ar_multiply_W_new';
            sigma_ky_new = zeros(L_new,1);
            for p = 1:Nx
                sigma_ky_new = sigma_ky_new + (Ar_multiply_W_new*Y(:,p)).*conj(At_multiply_X_new(:,p));
            end
            sigma_kk_new = At_multiply_X2_new.*Ar_multiply_W2_new;

            func_val = -sigma_ky_new'*((D/lambda + sigma_kk_new)\sigma_ky_new);
            maxit=maxit+1;
        end
        stepsize_old = stepsize;
        if (maxit < 30)
            theta = theta_new;
        end

        At = exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
        Ar = exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
        At_multiply_X = At'*X;
        At_multiply_X2 = (At_multiply_X*At_multiply_X').';
        Ar_multiply_W = Ar';
        Ar_multiply_W2 = Ar_multiply_W*Ar_multiply_W';
        sigma_ky = zeros(L_new,1);
        for p = 1:Nx
            sigma_ky = sigma_ky + (Ar_multiply_W*Y(:,p)).*conj(At_multiply_X(:,p));
        end
        sigma_kk = At_multiply_X2.*Ar_multiply_W2;
        
        z_new = (D/lambda + sigma_kk)\(sigma_ky);
    end
    z_old = z_new;
end
[z_es z_sort] = sort(z_old,'descend');
theta_es = theta(:,z_sort);
At = exp(-2i*pi*(0:1:Nt-1)'*theta(1,:));
Ar = exp(-2i*pi*(0:1:Nr-1)'*theta(2,:));
err = norm(Y-Ar*diag(z_es)*At'*X,'fro');
end