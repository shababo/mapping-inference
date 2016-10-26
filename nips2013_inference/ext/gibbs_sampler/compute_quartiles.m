function Qs=compute_quartiles(Sgs,mus,Cs)
% COMPUTE_QUARTILES computes the quartiles of the sample distribution of W.

nG=size(Sgs,2);
p=size(Sgs,1);
Qs=zeros(p,3);
qs=[1/4 1/2 3/4];
for q=1:length(qs) % for each quartile value
    for j=1:p % for each weight
        display(['Computing quartile ' num2str(q) ' for component ' num2str(j)])
        left=0; % \int_{-\infty}^0 p_{sample}(W|D)[Gaussians part]
        for g=1:nG % for each Gibbs sweep
            if Sgs(j,g)
                left=left+normcdf(0,mus(j,g),Cs(j,g));
            end
        end
        left=left/nG;
        middle=sum(Sgs(j,:)==0)/nG; % \int p_{sample}(W|D)[deltas part]
        if left<qs(q) && left+middle>qs(q)
            Qs(j,q)=0;
        else
            x=0;
            dx=100;
            err=1;
            while abs(err)>1e-3
                y=middle*nG*(x>0);
                for g=1:nG
                    if Sgs(j,g)
                        y=y+normcdf(x,mus(j,g),Cs(j,g));
                    end
                end
                y=y/nG;
                err=qs(q)-y;
                dx=sign(err)*abs(dx)/(1+(sign(err)~=sign(dx)));
                x=x+dx;
            end
            Qs(j,q)=x;
        end
    end
end
                     
end