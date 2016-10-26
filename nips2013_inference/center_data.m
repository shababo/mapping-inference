function [X_center, Y_center] = center_data(X, Y)
%CENTER_DATA(X, Y)
%
% pre-processing before varbvs: center the data in X and Y
%
% centering the Y vector: removes the intercept
%
% centering the X matrix: better approximation of sparse input
%


Y_bar = mean(Y);
Y_center = Y - Y_bar;
X_bar = mean(X,1);
X_center = X - X_bar(ones(length(Y), 1), :);