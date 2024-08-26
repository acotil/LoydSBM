Methods=["LloydForward","LloydBackward","LloydMix","LloydMLE","LloydTest","GD","VEM","Spectral","SpectralPy"]        # Methods is the vector of all admissible method to estimate a SBM from a LabeledGraph
Metrics=["l1","l2","Huber","linfty","Ent","Rent","l1c","l2c","linftyc","Entc","l1p","l2p","Huberp","linftyp","Entp","Rentp","l1cp","l2cp","linftycp","Entcp"]                 # Metrics is the vector of all admissible metric to use the LloydForward and the LloydBackward method in the estimation


"""
Structure representing the set of parameters for the estimation of an SBM from a LabeldGraph

Arguments
----------
method : String
    method used for the estimation.

with_ini : Bool
    If true then an initial labelisation must be given, if false then the initial labelisation is wompute from the spectral clusterinf method 

max_iter : Int
    maximum number of iteration if the stop condition is based on the number of iteration. Here one iteration correspond to the relabelisation of each node.

max_time : Float64
    maximum time allowed for the algorithm to calculate an SBM estimate.

metric : String
    metric used for the estimation if method is Lloyd type.

Huber_const : Union{Float64,Nothing}
    the parameter of the Huber loss function if Huber metric is used.
    
tol : Float64   
    tolerance threshold used in the stoping condition if the method is VEM

tol : Float64   
    tolerance threshold used in the stoping condition of the fixed point algorithm if the method is VEM

reg_coeff : Float64
    the coefficient of regularisation used if the method is Spectral
    
Vebose : Bool
    NOT YET IMPLEMENTED (must be set to false).
"""
struct SBMEstimator
    method::String                              # Method used
    
    # Parameters iterative methods
    with_ini::Union{Nothing,Bool}
    max_iter::Union{Nothing,Int}                # Max number of iter to converge
    max_time::Union{Nothing,Float64}            # Max time to converge
    
    # Parameters Lloyd
    metric::Union{Nothing,String}               # Parameters used for metric
    Huber_const::Union{Nothing,Float64}         # parameter of the Huber loss function

    # Parameters VEM
    tol::Union{Nothing,Float64}                 # Parameter used to define when the two lap_matrix are equal after one iteration of the VEM algorithm
    tolFP::Union{Nothing,Float64}               # Parameter used to define when the two lap_matrix are equal after one iteration of the Fixed-point algorithm
    
    # Parameters Spectral
    reg_coeff::Union{Nothing,Float64}           # regularisation parameter used in the spectral clustering
    
    # Optionnal general parameters 
    verbose::Bool                               # Verbose

    function SBMEstimator(method::String ;
                    with_ini::Union{Nothing,Bool} = false,
                    max_iter::Union{Nothing,Int} = nothing, 
                    max_time::Union{Nothing,Float64} = nothing, 
                    metric::Union{Nothing,String} = nothing, 
                    Huber_const::Union{Nothing,Float64} = nothing,
                    tol::Union{Nothing,Float64} = nothing,
                    tolFP::Union{Nothing,Float64} = nothing,  
                    reg_coeff::Union{Nothing,Float64} = nothing, 
                    verbose::Bool=false)::SBMEstimator
        
        if method in ["LloydForward","LloydBackward", "LloydMix","LloydMLE", "LloydTest", "GD"]

            if max_iter === nothing
                throw(ArgumentError("If method is $method, max_iter must be given"))
            elseif max_iter <= 0
                throw(ArgumentError("max_iter must be a positive value (in number of steps)"))
            end
    
            if max_time === nothing
                throw(ArgumentError("If method is $method, max_time must be given"))
            elseif max_time <= 0
                throw(ArgumentError("max_time must be a positive value (in number of seconds)"))
            end

            if method in ["LloydForward","LloydBackward", "LloydMix"]
                if metric === nothing
                    throw(ArgumentError("If method is $method, metric must be given"))
                elseif ~(metric in Metrics)
                    throw(ArgumentError("Invalid method: $metric."))
                end
    
                if metric == "Huber"
                    if Huber_const === nothing
                        throw(ArgumentError("If method is $method and metric is Huber, Huber_const must be given"))
                    elseif Huber_const <= 0.0
                        throw(ArgumentError("If method is LloydForward, LloydBackward or LloydMix and metric is Huber or Huberc, the Huber_const must be positive"))
                    end
                end
            end
        elseif method in ["VEM"]
            if max_iter === nothing
                throw(ArgumentError("If method is $method, max_iter must be given"))
            elseif max_iter <= 0
                throw(ArgumentError("max_iter must be a positive value (in number of steps)"))
            end
    
            if max_time === nothing
                throw(ArgumentError("If method is $method, max_time must be given"))
            elseif max_time <= 0
                throw(ArgumentError("max_time must be a positive value (in number of seconds)"))
            end

            if tol === nothing
                throw(ArgumentError("If method is $method, tol must be given"))
            elseif tol <= 0
                throw(ArgumentError("tol must be a positive value (in number of seconds)"))
            end

            if tolFP === nothing
                throw(ArgumentError("If method is $method, tolFP must be given"))
            elseif tolFP <= 0
                throw(ArgumentError("tolFP must be a positive value (in number of seconds)"))
            end
        elseif method in ["Spectral", "SpectralPy"]
            if reg_coeff !== nothing
                if reg_coeff <= 0
                    throw(ArgumentError("If method is Spectral or SpectralPy and if reg_coeff is given then it must be positive"))
                end
            end
        else
            throw(ArgumentError("Invalid method: $method."))
        end

        return new(method, with_ini, max_iter, max_time, metric, Huber_const, tol, tolFP, reg_coeff, verbose)
    end
end

include("SpectralOpt.jl")
include("GDOpt.jl")
include("LloydOpt.jl")
include("VEMOpt.jl")



"""
Estimate the labelisation from estimator.method of adj_matrix with a number of n_labels of labels

Parameters
----------
estimator : SBMEstimator
    A Estimator structure, giving the method and the related parameters of the estimation.

n_labels : Int
    The expected number of labels

adj_matrix : Matrix{Int}
    An adjacency matrix.

Optional Parameters
----------

ini_lab_vector : Union{Nothing,Vector{Int}}
    The initial vector of labels if estimator.with_ini is true

ini_lab_matrix : Union{Nothing,Matrix{Float64}}
    The initial matrix of labels if estimator.with_ini is true

seed : Union{Nothing,Int}
    Seed for the random number generator.

Returns
-------
new_lab_vector : Vector{Int}
    the estimate vector of labels

estim_lab_matrix : Matrix{Float64}
    the estimate matrix of labels 

estim_prob_labels : Vector{Float64}
    the estimate probability vector

 estim_prob_matrix : Matrix{Float64}
    the estimate matrix of probability 

estim_LLH : Float64
    the corresponding value of the log-likelyhood

computation_time : Float64
    the estimation time

n_iter : Int
    the number of iteration of Lloyd or VEM algorithm

status : String
    The status of the estimation, can  be CONVERGENCE if the algorithm converged, ITER_LIMIT if the maximum number of iterations has been reached, TIME_LIMIT if the maximum time has been reached or OTHER_STATUS if something went wrong.
"""
function Estimate(estimator::SBMEstimator, n_labels::Int, adj_matrix::Matrix{Int} ; ini_lab_vector::Union{Nothing,Vector{Int}} = nothing, ini_lab_matrix::Union{Nothing,Matrix{Float64}} = nothing, seed::Union{Nothing,Int}=nothing)
    
    # Set the random seed if one was provided
    if seed !== nothing
        Random.seed!(seed)
    end

    #----------------------------

    # Loading parameters

    method = estimator.method
    metric = estimator.metric
    max_iter = estimator.max_iter
    max_time = estimator.max_time
    Huber_const = estimator.Huber_const
    reg_coeff = estimator.reg_coeff
    tol = estimator.tol
    tolFP = estimator.tolFP
    with_ini = estimator.with_ini
    n_nodes = size(adj_matrix,1)

    # Loading the Optimizer related to the used method 

    if method in ["LloydForward", "LloydBackward", "LloydMix", "LloydMLE", "LloydTest", "GD"]
        
        if ~with_ini 
            ini_lab_vector = SpectralOpt(n_nodes, n_labels, adj_matrix, reg_coeff)[1]
            # Check if all labels are present
            if length(Set(ini_lab_vector)) != n_labels
                throw(ArgumentError("The initialisation found by the Septral method does not contain the right number of labels"))
            end
        else
            # Check if data is initialized
            if ini_lab_vector === nothing && ini_lab_matrix === nothing
                throw(ArgumentError("if with_ini = true then ini_lab_vector or ini_lab_matrix must be given"))
            elseif ini_lab_vector === nothing && ini_lab_matrix !== nothing
                ini_lab_vector = lab_matrix_to_vector(ini_lab_matrix)
            end
            # Check if all labels are present
            if ~(allpresent(ini_lab_vector,n_labels))
                throw(ArgumentError("All labels must be represented in ini_lab_vector (the initialization)."))
            end
        end
        
        return(LloydEstimate(n_labels, ini_lab_vector, adj_matrix, method, max_iter, max_time, metric, Huber_const))

    elseif method == "VEM"
        if ~with_ini 
            ini_lab_vector = SpectralOpt(n_nodes, n_labels, adj_matrix, reg_coeff)[1]
            ini_lab_matrix = lab_vector_to_matrix(n_labels, ini_lab_vector)
            # Check if all labels are present
            if length(Set(ini_lab_vector)) != n_labels
                throw(ArgumentError("The initialisation found by the Septral method does not contain the right number of labels"))
            end
        else
            # Check if data is initialized
            if ini_lab_vector === nothing && ini_lab_matrix === nothing
                throw(ArgumentError("if with_ini = true then ini_lab_vector or ini_lab_matrix must be given"))
            elseif ini_lab_matrix === nothing && ini_lab_vector !== nothing
                ini_lab_matrix = lab_vector_to_matrix(n_labels, ini_lab_vector)
            end
            # Check if all labels are present
            if ~(allpresent(ini_lab_matrix))
                throw(ArgumentError("All labels must be represented in ini_lab_matrix (the initialization)."))
            end
        end

        return(VEMEstimate(n_labels, ini_lab_matrix, adj_matrix, max_iter, max_time, tol, tolFP))

    elseif method == "Spectral"
        return(SpectralOpt(n_nodes, n_labels, adj_matrix, reg_coeff))
    elseif method == "SpectralPy"
        return(SpectralPyOpt(n_nodes, n_labels, adj_matrix, reg_coeff))
    else
        throw("no method matching with $method")
    end
end