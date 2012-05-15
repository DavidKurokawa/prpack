classdef pagerank_solver < handle
    %PAGERANK_SOLVER PageRank solver
    %   A PageRank solver that upon creation takes in a graph, and between
    %   subsequent calls of solve will store preprocessed work.
    
    properties (Hidden)
        solver_ptr;
    end
    
    methods
        % constructor
        function obj = pagerank_solver(matrix)
            obj.solver_ptr = create_pagerank_solver(matrix);
        end
        % destructor
        function delete(obj)
            delete_pagerank_solver(obj.solver_ptr);
        end
        % method
        function [x, ret] = solve(obj, alpha, tol, u, v, method)
            [x, ret] = use_pagerank_solver(obj.solver_ptr, alpha, tol, u, v, method);
        end
    end
    
end
