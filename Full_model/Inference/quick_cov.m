function [cov_out] = quick_cov(...
    A,B)
dA=A-mean(A);dB=B-mean(B);
cov_out=mean(dA.*dB);
end
