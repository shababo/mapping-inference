function L=full_conditional_sampler_lambda(b)

% f=@(x) 1./x*1./(1+x).*exp(-b./x); % target distribution, up to a constant
% g_1=@(x) 1./x.^2.*exp(-b./x); % one envelope
% g_2=@(x) 1./(2*x.^(3/2)).*exp(-b./x); % a complementary envelope

f=@(x) 1./x*1./(1+x).*exp(-b./x); % target distribution, up to a constant
g_1=@(x) 1./x.^2.*exp(-b./x); % one envelope
g_2=@(x) 1./(2*x.^(3/2)).*exp(-b./x); % a complementary envelope

watch=0; % Controls whether to plot rejection sampling in real time

if watch
    % plot envelopes against target function
    h_fig_1=figure(1);
    set(h_fig_1,'units','inches')
    set(h_fig_1,'outerposition',[0 2 15 10])
    clf
    subplot(2,1,1)
    x=linspace(0,10*b,100);
    hold on
    plot(x,f(x),'k')
    plot(x,g_1(x),'r')
    subplot(2,1,2)
    x=linspace(0,10*b,100);
    hold on
    plot(x,f(x),'k')
    plot(x,g_2(x),'r')
end

tries=0;
while 1
    tries=tries+1;
    % Rejection sample with g(x)=inversegamma(1,b)
    u=1/gamrnd(1,1/b);
    r=rand;
    if watch
        % Mark proposed L value vs r*f(u)/g(u)
        subplot(2,1,1)
        scatter(u,r*g_1(u),'x')
%         drawnow
    end
    if r<f(u)/g_1(u)
        L=u;
        break;
    end
    % Rejection sample with g(x)=2*inversegamma(1/2,b)
    u=1/gamrnd(1/2,1/b);
    r=rand;
    if watch
        % Mark proposed L value vs r*f(u)/g(u)
        subplot(2,1,2)
        scatter(u,r*g_2(u),'x')
%         drawnow
    end
    if r<f(u)/g_2(u)
        L=u;
        break;
    end
end
% disp([num2str(tries) ' tries to sample lambda, b = ' num2str(b)])