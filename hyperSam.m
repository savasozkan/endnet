function [errRadians] = hyperSam(a, b)
% HYPERSAM Computes the spectral angle error (in radians) between two
% vectors, or between a vector and every column or row of a matrix
%
% Usage
%   [errRadians] = hyperSam(a, b)
%   [errRadians] = hyperSam(A, b)
%   [errRadians] = hyperSam(a, B)
% Inputs
%   a - vector, (px1)
%   b - vector, (px1)
%   A - matrix, (pxN)
%   B - matrix, (pxN)
% Outputs
%   errRadians - angle between vectors a and b in radians, (1xN)

% Check dimensions

% Turn row vectors to column vectors if necessary

    errRadians = getSam(a,b);


function errRadians = getSam(a, b)
% GETSAM computes the spectral angle error between the vector b and
% every column of matrix a (or just the first column, if a is a vector)

[~,N] = size(a);
[~,L] = size(b);
errRadians = zeros(L,N);

for l=1:L
    for k=1:N
        tmpa = a(:,k);
        tmpb = b(:,l);
        errRadians(l,k) = acos( dot(tmpa, tmpb)/(norm(tmpb)*norm(tmpa)) );
    end
end