function samples = rnd_truncated_normal(truncs)

% Author: Ari Pakman, using technique from Robert's 2009 paper.
% modified by: Ben Shababo 2013.05.10

% returns samples from a normal distribution with zero mean and unit
% variance, truncated to be bigger than truncs

% Input 
% truncs = vector of size (N,1) with the points of truncation. 

% Output
% sample:   a vector with the samples 



N = max(size(truncs));
samples = zeros(N,1);

% for negative truncation, sample from untruncated normal until the
% result falls to the right of the truncation
tosample = (truncs<=0);
tsn = sum(tosample);

while tsn
    zs=normrnd(0,1, tsn,1);
    
    keep = false(N,1);
    tokeep= (zs >= truncs(tosample));
    
    keep(tosample) = tokeep;
    samples(keep)=zs(tokeep);
    tosample(keep)=false;
    tsn = sum(tosample);
end %while




% for positive truncation, use accept-reject with a an exponential
% distribution as bound

tosample = (truncs > 0);

alphas = 0.5*(truncs + sqrt(truncs.^2 + 4));
tsp = sum(tosample);

while tsp
   tsp;
    zs=random('exp',1./alphas(tosample)) + truncs(tosample);
    rhos = exp(-0.5*(zs-alphas(tosample)).^2);
    us = rand(1,tsp);

    keep = false(N,1);
    tokeep= (us<=rhos);
    keep(tosample) = tokeep;
    samples(keep)=zs(tokeep);
    tosample(keep)=false;
    tsp = sum(tosample);
end %while

end

