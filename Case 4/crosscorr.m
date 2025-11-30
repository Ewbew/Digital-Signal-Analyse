
function [r, lags] = crosscorr(x, y)
% SIMPLE_XCORR  Minimal cross-correlation of two vectors
%   [r, lags] = SIMPLE_XCORR(x, y) returns r(k) for lags from -(Ny-1) to (Nx-1)
%   where Nx = length(x), Ny = length(y). Works for real or complex vectors.
%
%   Example:
%     x = [1 2 3];
%     y = [4 5];
%     [r,lags] = simple_xcorr(x,y);

x = x(:).';    % ensure row vectors
y = y(:).';

% flip y and convolve
r_full = conv(x, fliplr(y));

% lags: -(length(y)-1) : (length(x)-1)
lags = -(length(y)-1) : (length(x)-1);

r = r_full;    % already in correct order for these lags
end