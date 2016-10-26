function [X_s, w_s, gamma_s] = gibbs_single_sweep_known_delays(X_s, D, w_s, gamma_s, Y, pi_nk, c, a, sigma_s, sigma_n, d_mean_nk, d_sigma_nk, t, tau, gmax)

[N, K] = size(X_s);

% let's do X and D
fprintf('sampling X and D\n');
tic
%[X_s, D_s] = slow_sample_X_D(Y,X_s,D_s,w_s,pi_nk,d_mean_nk,d_sigma_nk,t,tau,gmax,sigma_n);
[X_s, D_s] = sample_X(Y,double(X_s),D,w_s,pi_nk,d_mean_nk,d_sigma_nk,t,tau,gmax,sigma_n);
X_s = X_s > 0;
toc

% compute sufficient statistics M and r
fprintf('computing sufficient statistics\n');
tic
[M, r] = compute_M_and_r(Y, X_s, D, t, tau, gmax, sigma_n, ones(K,1));
toc

% test if M is symmetric and pos-semidef
if any(eig(M) < 0)
    fprintf('M is not pos semi-def :( \n')
end

M  = M + diag(1./(sigma_s.^2));
% T = size(Y,2);
% r_test = zeros(K,1);
% M_test = zeros(K);
% D_s_t = zeros(N,length(t),K);
% for k = 1:K
% 	D_s_t(:,:,k) = alpha_synapse(t,D_s(:,k),tau,-gmax);
% 	r_test(k) = sum(X_s(:,k)'*(D_s_t(:,:,k).*Y))/sigma_n^2;
% 	M_test(k,k) = sum(sum((X_s(:,k*ones(T,1)).*D_s_t(:,:,k)).^2))/sigma_n^2;
% end
% 
% for k = 1:K
%     for j = k+1:K
%         M_test(k,j) = sum(sum((X_s(:,k*ones(T,1)).*D_s_t(:,:,k)).*(X_s(:,j*ones(T,1)).*D_s_t(:,:,j))))/sigma_n^2;
%         M_test(j,k) = M_test(k,j);
%     end
% end
% 
% disp('M diff')
% sum(norm(M(:) - M_test(:)))
% disp('r diff')
% sum(norm(r(:) - r_test(:)))



fprintf('sampling w and gamma\n');
tic
% sample w and \gamma (Ari Pakman)
for k = 1:K

    gamma_s(k) = 0;
    mask = gamma_s > 0;      
    mi = M(k,mask);

    alpha = M(k,k);
    beta = r(k) - mi*w_s(mask);

    q = sign(c(k)-.5)*beta/sqrt(alpha);

    ncdfq= normcdf(q,0,1);

    if normcdf(q,0,1) == 0
        odds_ratio = 0;
    else
        odds_ratio = exp(q^2/2) * (a(k)/(1-a(k))) * 2/(sigma_s(k)*sqrt(alpha)) * ncdfq;
    end

    if odds_ratio == Inf || rand() < odds_ratio/(odds_ratio+1)
        gamma_s(k) = 1;
        w_k = rnd_truncated_normal(-q);
        w_s(k) = sign(c(k)-.5)*(w_k + q)/sqrt(alpha);
    else
        gamma_s(k) = 0;                                
        w_s(k) = 0;

    end
    
    if isinf(w_s(k))
        disp('doh, a w sample is Inf!')
        return
    end

end
toc
