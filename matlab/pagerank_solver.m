classdef pagerank_solver < handle
    %PAGERANK_SOLVER PageRank solver
    %   A PageRank solver that upon creation takes in a graph, and between
    %   subsequent calls of solve will store preprocessed work.
    
    properties (Hidden)
        solver;
    end
    
    methods
        function obj = pagerank_solver(num_vs, heads, tails)
            obj.solver = create_pagerank_solver(num_vs, heads, tails);
        end
        function [x, ret] = solve(obj, alpha, tol, u, v, method)
            [x, ret] = use_pagerank_solver(obj.solver, alpha, tol, u, v, method);
        end
    end
    
end

