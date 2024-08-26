"""
Given a labelisation of nodes, compute the relabelisation that minimizes the LLH by relabelised once each node and keeping the labels of the other nodes unchanged

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
function GDOpt(n_nodes::Int, n_labels::Int, lab_vector::Vector{Int}, adj_matrix::Matrix{Int}, metric::Union{String,Nothing}=nothing, Huber_const::Union{Float64,Nothing}=nothing)

    estim_lab_vector = zeros(Int,n_nodes)
    
    for i in 1:n_nodes
        current_label = lab_vector[i]
        estim_LLH = - Inf
        for k in 1:n_labels
            lab_vector[i] = k
            prob_matrix = CalculateProbMatrix(n_labels, lab_vector, adj_matrix)
            LLH = CalculateLLH(lab_vector, prob_matrix, adj_matrix)
            if LLH > estim_LLH
                estim_LLH = LLH
                estim_lab_vector[i] = k 
            end
        end
        lab_vector[i] = current_label
    end

    return(estim_lab_vector)

end