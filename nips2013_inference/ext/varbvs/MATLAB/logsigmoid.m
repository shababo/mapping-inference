% LOGSIGMOID(X) is the logarithm of the sigmoid of the elements of X. Use
% this instead of LOG(SIGMOID(X)) to avoid loss of numerical precision.
function y = logsigmoid (x)
  y = -logpexp(-x);
  