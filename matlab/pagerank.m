function [x, stats] = pagerank(matrix, alpha, tol, u, v, method)
%PAGERANK Computes the PageRank
%   Computes the PageRank of the sparse matrix given with the given
%   parameters.
%   Inputs:
%       matrix = a sparse matrix whose non-zero entries indicate edges (all
%                entry values are ignored)
%       alpha  = alpha value
%       tol    = tolerance allowed for error
%       u      = u vector
%       v      = v vector
%       method = [optional] method to use to compute PageRank

prs = pagerank_solver(matrix);
if (nargin < 4)
    error('Too few input arguments.');
elseif (nargin == 5)
    [x, stats] = prs.solve(alpha, tol, u, v);
elseif (nargin == 6)
    [x, stats] = prs.solve(alpha, tol, u, v, method);
else
    error('Too many input arguments.');

end
