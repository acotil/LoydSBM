"""
Given a labelisation of nodes, compute the relabelisation that minimizes the intra-class variance (related to the given metric) of the empirical emission laws by keeping the center of each class

Parameters
----------
metric : String
    the metric used to compute the intra-class variance
n_nodes : Int
    the number of nodes of the LabeledGraph
n_labels : Int
    the expected nuber of labels
lab_vector : Vector{Int}
    the initial vector of labels
adj_matrix : Matrix{Int}
    the adjacency matrix of the Labeled GraphLabeled
Huber_const : Float64
    the parameter of the Huber loss function if metric is Huber

Returns
-------
estim_lab_vector : Vector{Int}
    The new labelisation
estim_prob_matrix : Matrix{Float64}
    The new probability matrix computed from the new labelisation
estim_LLH : Float64
    The value of the log-likelyhood computed from the new labelisation and the related probability matrix
"""
function LloydForwardOpt(n_nodes::Int, n_labels::Int, lab_vector::Vector{Int}, adj_matrix::Matrix{Int}, metric::String, Huber_const::Union{Float64,Nothing} = nothing)
        
    Z=[Int(l == k) for l in lab_vector, k in 1:n_labels]
    N = vec(sum(Z, dims=1))
        
    # Computation of estim_pi[i,k], the estimated probability that a node i is linked to a node of label k 
    x = adj_matrix * Z                                  # x[i,l] is the number of edges from i to nodes of label l
    estim_pi = zeros(n_nodes,n_labels)                  # estim_pi[i,l] is the estimated probability that a node i is linked to a node of label l
    for i in 1:n_nodes, l in 1:n_labels
        if (N[l] != 0)                              # If there is no nodes of label l then estim_pi[i,l]=0 for all i
            estim_pi[i, l] = x[i, l] / N[l]
        end
    end

    # Computation of estim_prob_matrix[k,l], the estimated probability that a node of label k is linked to a node of label l
    y = transpose(Z) * estim_pi
    estim_prob_matrix = zeros(n_labels,n_labels)                  # estim_prob_matrix[k,l] is the estimated probability that a node of label k is linked to a node of label l
    for k in 1:n_labels, l in 1:n_labels
        if (N[k] != 0)                              # If there is no nodes of label k then estim_prob_matrix[k,l]=0 for all l
            estim_prob_matrix[k, l] = y[k, l] / N[k]
        end
    end

    #println(CalculateObjective(metric, n_labels, n_nodes, labels, estim_pi, estim_prob_matrix, x, N, Huber_const, pen))

    estim_lab_vector = [argmin([Dist(metric, estim_pi[i, :], estim_prob_matrix[k, :], x[i, :], N, Huber_const,Float64(N[k]/n_labels)) for k in 1:n_labels]) for i in 1:n_nodes]

    return(estim_lab_vector)
end

"""
Given a labelisation of nodes, compute the relabelisation that minimizes the intra-class variance (related to the given metric) of the empirical reception laws by keeping the center of each class

Parameters
----------
metric : String
    the metric used to compute the intra-class variance
n_nodes : Int
    the number of nodes of the LabeledGraph
n_labels : Int
    the expected nuber of labels
lab_vector : Vector{Int}
    the initial vector of labels
adj_matrix : Matrix{Int}
    the adjacency matrix of the Labeled GraphLabeled
Huber_const : Float64
    the parameter of the Huber loss function if metric is Huber

Returns
-------
estim_lab_vector : Vector{Int}
    The new labelisation
estim_prob_matrix : Matrix{Float64}
    The new probability matrix computed from the new labelisation
estim_LLH : Float64
    The value of the log-likelyhood computed from the new labelisation and the related probability matrix
"""
function LloydBackwardOpt(n_nodes::Int, n_labels::Int, lab_vector::Vector{Int}, adj_matrix::Matrix{Int}, metric::String, Huber_const::Union{Float64,Nothing} = nothing)
        
    Z=[Int(l == k) for l in lab_vector, k in 1:n_labels]
    N = vec(sum(Z, dims=1))
    
    # Computation of estim_pi[i,k], the estimated probability that a node i is linked to a node of label k 
    x = transpose(Z) * adj_matrix                                # x[i,l] is the number of edges from i to nodes of label l
    estim_pi = zeros(n_labels,n_nodes)                  # estim_pi[i,l] is the estimated probability that a node i is linked to a node of label l
    for i in 1:n_nodes, k in 1:n_labels
        if (N[k] != 0)                              # If there is no nodes of label l then estim_pi[i,l]=0 for all i
            estim_pi[k, i] = x[k, i] / N[k]
        end
    end

    # Computation of estim_prob_matrix[k,l], the estimated probability that a node of label k is linked to a node of label l
    y = estim_pi * Z
    estim_prob_matrix = zeros(n_labels,n_labels)                  # estim_prob_matrix[k,l] is the estimated probability that a node of label k is linked to a node of label l
    for k in 1:n_labels, l in 1:n_labels
        if (N[l] != 0)                              # If there is no nodes of label k then estim_prob_matrix[k,l]=0 for all l
            estim_prob_matrix[k, l] = y[k, l] / N[l]
        end
    end

    #println(CalculateObjective(metric, n_labels, n_nodes, lab_vector, transpose(estim_pi), transpose(estim_prob_matrix), y, N, Huber_const, pen))

    estim_lab_vector = [argmin([Dist(metric, estim_pi[:, i], estim_prob_matrix[:, k], x[:, i], N, Huber_const,Float64(N[k]/n_labels)) for k in 1:n_labels]) for i in 1:n_nodes]

    return(estim_lab_vector)
end

"""
Given a labelisation of nodes, compute the relabelisation that minimizes the intra-class variance (related to the given metric) of the empirical emission and reception laws by keeping the center of each class

Parameters
----------
metric : String
    the metric used to compute the intra-class variance
n_nodes : Int
    the number of nodes of the LabeledGraph
n_labels : Int
    the expected nuber of labels
lab_vector : Vector{Int}
    the initial vector of labels
adj_matrix : Matrix{Int}
    the adjacency matrix of the Labeled GraphLabeled
Huber_const : Float64
    the parameter of the Huber loss function if metric is Huber

Returns
-------
estim_lab_vector : Vector{Int}
    The new labelisation
estim_prob_matrix : Matrix{Float64}
    The new probability matrix computed from the new labelisation
estim_LLH : Float64
    The value of the log-likelyhood computed from the new labelisation and the related probability matrix
"""
function LloydMixOpt(n_nodes::Int, n_labels::Int, lab_vector::Vector{Int}, adj_matrix::Matrix{Int}, metric::String, Huber_const::Union{Float64,Nothing} = nothing)
        
    Z=[Int(l == k) for l in lab_vector, k in 1:n_labels]
    N = vec(sum(Z, dims=1))
    
    # Computation of estim_pi[i,k], the estimated probability that a node i is linked to a node of label k 
    x = adj_matrix * Z                                  # x[i,l] is the number of edges from i to nodes of label l
    estim_out = zeros(n_nodes,n_labels)                  # estim_pi[i,l] is the estimated probability that a node i is linked to a node of label l
    for l in 1:n_labels
        if (N[l] != 0)                              # If there is no nodes of label l then estim_pi[i,l]=0 for all i
            estim_out[:, l] = x[:, l] / N[l]
        end
    end

    # Computation of estim_nu[k,j], the estimated probability that a node of label k is linked to a node j 
    y = transpose(Z) * adj_matrix                                # x[i,l] is the number of edges from i to nodes of label l
    w = transpose(Z) * estim_out
    estim_in = zeros(n_labels,n_nodes)                           # estim_pi[i,l] is the estimated probability that a node i is linked to a node of label l
    estim_prob_matrix = zeros(n_labels,n_labels)                  # estim_prob_matrix[k,l] is the estimated probability that a node of label k is linked to a node of label l
    for k in 1:n_labels
        if (N[k] != 0)  
            estim_in[k, :] = y[k, :] / N[k]
            estim_prob_matrix[k, :] = w[k, :] / N[k]
        end
    end

    #estim_lab_vector = [argmin([Dist(metric, vcat(estim_out[i, :], estim_in[:,i]), vcat(estim_prob_matrix[k, :], estim_prob_matrix[:, k]), vcat(x[i, :], y[:, i]), N, Huber_const,Float64(N[k]/n_labels)) for k in 1:n_labels]) for i in 1:n_nodes]

    estim_lab_vector = zeros(Int,n_nodes)
    for i in 1:n_nodes
        estim_label = 0
        estim_obj = Inf
        for k in 1:n_labels
            obj = Dist(metric, vcat(estim_out[i, :], estim_in[:,i]), vcat(estim_prob_matrix[k, :], estim_prob_matrix[:, k]), vcat(x[i, :], y[:, i]), N, Huber_const,Float64(N[k]/n_labels))
            if obj < estim_obj
                estim_obj = obj
                estim_label = k
            end
        end
        estim_lab_vector[i] = estim_label
    end
    return(estim_lab_vector)
end

"""
Given a labelisation of nodes, compute the relabelisation by method for test 

Parameters
----------
metric : String
    the metric used to compute the intra-class variance
n_nodes : Int
    the number of nodes of the LabeledGraph
n_labels : Int
    the expected nuber of labels
lab_vector : Vector{Int}
    the initial vector of labels
adj_matrix : Matrix{Int}
    the adjacency matrix of the Labeled GraphLabeled
Huber_const : Float64
    the parameter of the Huber loss function if metric is Huber

Returns
-------
estim_lab_vector : Vector{Int}
    The new labelisation
estim_prob_matrix : Matrix{Float64}
    The new probability matrix computed from the new labelisation
estim_LLH : Float64
    The value of the log-likelyhood computed from the new labelisation and the related probability matrix
"""
function LloydMLEOpt(n_nodes::Int, n_labels::Int, lab_vector::Vector{Int}, adj_matrix::Matrix{Int}, metric::Union{String,Nothing}=nothing, Huber_const::Union{Float64,Nothing}=nothing)
        
    estim_prob_matrix = CalculateProbMatrix(n_labels, lab_vector, adj_matrix)
    estim_lab_vector = zeros(Int,n_nodes)
    
    for i in 1:n_nodes
        estim_obj =  -Inf
        estim_label = lab_vector[i]
        for k in 1:n_labels
            obj = 0.0
            for j in 1:n_nodes
                if j!=i
                    if (adj_matrix[i,j]==1)
                        obj += log(estim_prob_matrix[k,lab_vector[j]])
                    else
                        obj += log(1-estim_prob_matrix[k,lab_vector[j]])
                    end
                    if (adj_matrix[j,i]==1)
                        obj += log(estim_prob_matrix[lab_vector[j],k])
                    else
                        obj += log(1-estim_prob_matrix[lab_vector[j],k])
                    end
                else 
                    if (adj_matrix[i,i]==1)
                        obj += 2*log(estim_prob_matrix[k,k])
                    else
                        obj += 2*log(1-estim_prob_matrix[k,k])
                    end
                end
            end
            if obj > estim_obj
                estim_obj = obj
                estim_label = k
            end
        end
        estim_lab_vector[i] = estim_label
    end

    return(estim_lab_vector)
end


"""
Given a labelisation of nodes, compute the relabelisation by method for test 

Parameters
----------
metric : String
    the metric used to compute the intra-class variance
n_nodes : Int
    the number of nodes of the LabeledGraph
n_labels : Int
    the expected nuber of labels
lab_vector : Vector{Int}
    the initial vector of labels
adj_matrix : Matrix{Int}
    the adjacency matrix of the Labeled GraphLabeled
Huber_const : Float64
    the parameter of the Huber loss function if metric is Huber

Returns
-------
estim_lab_vector : Vector{Int}
    The new labelisation
estim_prob_matrix : Matrix{Float64}
    The new probability matrix computed from the new labelisation
estim_LLH : Float64
    The value of the log-likelyhood computed from the new labelisation and the related probability matrix
"""
function LloydTestOpt(n_nodes::Int, n_labels::Int, lab_vector::Vector{Int}, adj_matrix::Matrix{Int}, metric::Union{String,Nothing}=nothing, Huber_const::Union{Float64,Nothing}=nothing)
    
    estim_lab_vector = copy(lab_vector)
    
    for i in 1:n_nodes
        estim_obj =  - Inf
        estim_label = lab_vector[i]
        estim_prob_matrix = CalculateProbMatrix(n_labels, estim_lab_vector, adj_matrix)
        
        for k in 1:n_labels
            obj = 0.0
            for j in 1:n_nodes
                if j!=i
                    if (adj_matrix[i,j]==1)
                        obj += log(estim_prob_matrix[k,estim_lab_vector[j]])
                    else
                        obj += log(1-estim_prob_matrix[k,estim_lab_vector[j]])
                    end
                    if (adj_matrix[j,i]==1)
                        obj += log(estim_prob_matrix[estim_lab_vector[j],k])
                    else
                        obj += log(1-estim_prob_matrix[estim_lab_vector[j],k])
                    end
                else 
                    if (adj_matrix[i,i]==1)
                        obj += 2*log(estim_prob_matrix[k,k])
                    else
                        obj += 2*log(1-estim_prob_matrix[k,k])
                    end
                end
            end
            if obj > estim_obj
                estim_obj = obj
                estim_label = k
            end
        end
        estim_lab_vector[i] = estim_label
    end
    
    return(estim_lab_vector)       
    
end


function LloydEstimate(n_labels::Int, ini_lab_vector::Vector{Int}, adj_matrix::Matrix{Int}, method::String, max_iter::Int, max_time::Float64, metric::Union{Nothing,String}=nothing, Huber_const::Union{Nothing,Float64}=nothing)
    
    if method == "LloydForward"
        LocalOpt = LloydForwardOpt
    elseif method == "LloydBackward"
        LocalOpt = LloydBackwardOpt
    elseif method == "LloydMix"
        LocalOpt = LloydMixOpt
    elseif method == "LloydMLE"
        LocalOpt = LloydMLEOpt
    elseif method == "LloydTest"
        LocalOpt = LloydTestOpt
    elseif method == "GD"
        LocalOpt = GDOpt
    else
        throw("no method matching with $method")
    end

    # Initialisation

    n_nodes = size(adj_matrix,1)
    n_iter = 0
    status = "OTHER_STATUS"
    start = time()
 
    # Heredity

    current_lab_vector = copy(ini_lab_vector)
    new_lab_vector = Vector{Int}()

    while true

        # Find the best relabelisation according to the used method        
        new_lab_vector = LocalOpt(n_nodes, n_labels, current_lab_vector, adj_matrix, metric, Huber_const)

        if islabequal(n_labels, current_lab_vector, new_lab_vector)
            status = "CONVERGENCE" 
            break
        elseif n_iter >= max_iter
            status = "ITER_LIMIT"
            break
        elseif time() - start >= max_time
            status = "TIME_LIMIT"
            break
        else
            n_iter += 1
            current_lab_vector = copy(new_lab_vector)
        end
    end
    
    estim_lab_matrix = lab_vector_to_matrix(n_labels, new_lab_vector)
    estim_prob_labels = CalculateProbLabels(n_labels, new_lab_vector)
    estim_prob_matrix = CalculateProbMatrix(n_labels, new_lab_vector, adj_matrix)
    estim_LLH = CalculateLLH(new_lab_vector, estim_prob_matrix, adj_matrix)
    
    computation_time = time()-start
    
    return(ini_lab_vector,new_lab_vector, estim_lab_matrix, estim_prob_labels, estim_prob_matrix, estim_LLH, computation_time, n_iter, status)
end

 