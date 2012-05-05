classdef pagerank_solver
    %PAGERANK_SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        solver;
    end
    
    methods
        function obj = pagerank_solver(num_vs, heads, tails)
            obj.solver = create_pagerank_solver(num_vs, heads, tails);
        end
        function ret = solve(obj, alpha, tol, u, v, method)
            [x, ret, obj.solver] = use_pagerank_solver(obj.solver, alpha, tol, u, v, method);
        end
    end
    
end

