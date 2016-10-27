function L=full_conditional_sampler_sigma(a,b)

f=@(x) x.^a*1./(1+x).*exp(-b*x); % target distribution, up to a constant
g_1=@(x) b^a/gamma(a)*x.^(a-1).*exp(-b*x); % one envelope
g_2=@(x) b^(a+1)/gamma(a+1)*x.^a.*exp(-b*x); % a complementary envelope

watch=0; % Controls whether to plot rejection sampling in real time

if watch
    % plot envelopes against target function
    h_fig_1=figure(1);
    clf
    x=linspace(0,10*b,100);
    subplot(1,2,1)
    hold on
    title('$\sigma$','interpreter','latex')
    plot(x,f(x),'k')
    plot(x,g_1(x)/(b^a/gamma(a)),'r')
    subplot(1,2,2)
    hold on
    plot(x,f(x),'k')
    plot(x,g_2(x)/(b^(a+1)/gamma(a+1)),'r')
end

tries=0;
while 1
    tries=tries+1;
    % Rejection sample with g(x)=gamma(a,b)
    u=gamrnd(a,1/b);
    f_over_g=f(u)*b^a/gamma(a)/g_1(u);
    r=rand;
    if watch
        % Mark proposed L value vs r*f(u)/g(u)
        figure(h_fig_1)
        subplot(1,2,1)
        scatter(u,r*g_1(u)*f_over_g,'x')
    end
    if r<f_over_g
        L=u;
        break;
    end
    % Rejection sample with g(x)=gamma(a-1,b)
    u=gamrnd(a+1,1/b);
    f_over_g=f(u)*(b^(a+1)/gamma(a+1))/g_2(u);
    r=rand;
    if watch
        % Mark proposed L value vs r*f(u)/g(u)
        figure(h_fig_1)
        subplot(1,2,2)
        scatter(u,r*g_2(u)*f_over_g,'x')
    end
    if r<f_over_g
        L=u;
        break;
    end
end
% disp([num2str(tries) ' tries to sample sigma'])