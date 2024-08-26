
"""
UNFINISHED
Estimate the labelisation from estimator.method of data
Parameters
----------
estimator : Estimator
    A Estimator structure, giving the method and the related parameters of the estimation.
data : LabeledGraph
    A LabeledGraph structure, giving the adjacency matrix and the initial labelisation.
seed : Union{Nothing,Int}
    Seed for the random number generator.

Returns
-------
estim_labels : Vector{Int}
    An estimation of the labelisation
estim_prob_matrix : Matrix{Float64}
    An estimation of prob_matrix
estim_LLH : Float64
    the value of the log-likelihood calculated from estim_labels, estim_prob_matrix and adj_matrix
computation_time : Float64
    the time that the computation to took
"""
function MINLP(data::LabeledGraph, time_limit::Float64=600.0, symm_break_const::Bool = false, symm_w_const::Bool =  false, w_lower::Float64 = 0.0)

    A = data.adj_matrix
    K = data.n_labels
    n = data.n_nodes

    # Upper bound on w
    w_upper = 1.0 - w_lower

    #Set time limit
    #= open("couenne.opt", "w") do f
        write(f, "time_limit $time_limit\n");
        write(f, "allowable_gap 1e-6\n");
        write(f, "allowable_fraction_gap 1e-6\n");
        write(f, "acceptable_tol 1e-6\n");
        write(f, "tol 1e-6\n");
        write(f, "feas_tolerance 1e-7\n");
        write(f, "integer_tolerance 1e-7\n");
        write(f, "constr_viol_tol 1e-7\n");
        write(f, "compl_inf_tol 1e-7\n"); 
        write(f, "acceptable_constr_viol_tol 1e-7\n");
        write(f, "acceptable_compl_inf_tol 1e-7\n");
    end =#

    # --- Optimization model -------------------------------------------------
    model = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe))
    

    # Variables
    @variable(model, w[1:K,1:K])
    @variable(model, z[1:n,1:K], Bin)

    # Objective
    @NLobjective(model, Max, 
           sum((log(w[r,s]) * z[i,r] * z[j,s]) for r = 1:K, s = 1:K, i = 1:n, j = 1:n if A[i,j] == 1) 
        + sum((log(1-w[r,s])* z[i,r] * z[j,s]) for r = 1:K, s = 1:K, i = 1:n, j = 1:n if A[i,j] == 0)
        )

    # Constraints
    @constraint(model, conwl[k = 1:K, l = 1:K], w[k,l] >= w_lower)
    @constraint(model, conwu[k = 1:K, l = 1:K], w[k,l] <= w_upper)
    @constraint(model, conz[i = 1:n], sum(z[i,r] for r = 1:K) == 1) # assignment constraint

    # Symmetry breaking in numbered clustering - Frank Plastria
    if symm_break_const
        # object 1 must be in cluster 1
        @constraint(model, z[1,1] == 1)
        @constraint(model, symBreak[r=2:(K-1), j=r:n], sum(z[i,l] for l = 1:(r-1), i = 2:(j-1)) - sum(z[j,l] for l = 1:r) <= j - 3 )
    end

    # Enforce symmetry on w
    if symm_w_const
        @constraint(model, symL[r=1:K, s=1:(r-1)], w[r,s] - w[s,r] <= 0.0)
        @constraint(model, symG[r=1:K, s=1:(r-1)], w[r,s] - w[s,r] >= 0.0)
    end
    # ------------------------------------------------------------------------

    # --- Solve optimization model -------------------------------------------
    
    optimize!(model)
    status = string(termination_status(model))
    println(status)
    
    #= if status == "LOCALLY_SOLVED"
        status = "OPTIMAL"
    elseif status == "OTHER_LIMIT"
        status = "TIME_LIMIT"
    end =#

    # ------------------------------------------------------------------------

    # --- Recover variable values --------------------------------------------
    @show w
    w_opt = value.(w)
    z_opt = value.(z)
    z_opt = round.(z_opt[:,:])
    z_opt[map(v -> isnan(v) , z_opt)] .= 0 # replace NaN values with 0
    z_opt = Int.(z_opt)
    # ------------------------------------------------------------------------

    println(z_opt)
    println(w_opt)

    return(z_opt, w_opt)
end
