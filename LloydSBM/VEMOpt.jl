
function localVEMOpt(n_nodes::Int, n_labels::Int, lab_matrix::Matrix{Float64}, adj_matrix::Matrix{Int}, tol_FP::Float64)
    
    # Computation of prob_labels at the n+1 step
    current_prob_labels = CalculateProbLabels(lab_matrix)

    # Computation of prob_matrix at the n+1 step
    current_prob_matrix = CalculateProbMatrix(lab_matrix,adj_matrix)

    # Computation of prob_apost at the n+1 step
    current_lab_matrix = reduce(vcat,transpose.([rand(Dirichlet(n_labels,1.0)) for i in 1:n_nodes]))         # ATTENTION : PAS CLAIR QUE CA SOIT CA QU'IL FAILLE FAIRE (essayer avec lab_matrix ou lab_matrix adoussie)
    new_lab_matrix = Matrix{Float64}([;;])
    n_iter = 0
    while true
        new_lab_matrix = ones(Float64, n_nodes, n_labels)
        n_iter += 1
        for i in 1:n_nodes
            for k in 1:n_labels
                for j in 1:n_nodes, l in 1:n_labels
                    a = 1.0
                    if adj_matrix[i,j] == 1
                        a *= current_prob_matrix[k,l]
                    else
                        a *= 1.0 - current_prob_matrix[k,l]
                    end
                    if adj_matrix[j,i] == 1
                        a *= current_prob_matrix[l,k]
                    else
                        a *= 1.0 - current_prob_matrix[l,k]
                    end
                    new_lab_matrix[i,k] *= a^(current_lab_matrix[j,l])
                end
                new_lab_matrix[i,k] *= current_prob_labels[k]
            end
            new_lab_matrix[i,:]./=sum(new_lab_matrix[i,:])
        end
        if islabequal(current_lab_matrix, new_lab_matrix, tol_FP) || n_iter == 10
            return(new_lab_matrix)
        else
            current_lab_matrix = copy(new_lab_matrix)
        end
    end
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
labels : Vector{Int}
    the initial vector of labels
adj_matrix : Matrix{Int}
    the adjacency matrix of the Labeled GraphLabeled
Huber_const : Float64
    the parameter of the Huber loss function if metric is Huber

Returns
-------
estim_labels : Vector{Int}
    The new labelisation
estim_prob_matrix : Matrix{Float64}
    The new probability matrix computed from the new labelisation
estim_LLH : Float64
    The value of the log-likelyhood computed from the new labelisation and the related probability matrix
"""
function VEMEstimate(n_labels::Int, ini_lab_matrix::Matrix{Float64}, adj_matrix::Matrix{Int}, max_iter::Int, max_time::Float64, tol::Float64, tolFP::Float64)

    # Initialisation

    n_iter = 0
    status = "OTHER_STATUS"
    start = time()
    n_nodes = size(ini_lab_matrix,1)

    # Heredity

    current_lab_matrix = copy(ini_lab_matrix)
    new_lab_matrix = Matrix{Float64}([;;])

    while true

        # Find the best relabelisation according to the used method        
        new_lab_matrix  = localVEMOpt(n_nodes, n_labels, current_lab_matrix, adj_matrix, tolFP)

        if islabequal(current_lab_matrix, new_lab_matrix, tol)
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
            current_lab_matrix = copy(new_lab_matrix)
        end
    end

    estim_lab_vector= lab_matrix_to_vector(new_lab_matrix)
    estim_prob_labels = CalculateProbLabels(new_lab_matrix)
    estim_prob_matrix = CalculateProbMatrix(new_lab_matrix, adj_matrix)
    estim_LLH = CalculateLLH(estim_lab_vector, estim_prob_matrix, adj_matrix)
    
    computation_time = time()-start
    
    return(lab_matrix_to_vector(ini_lab_matrix), estim_lab_vector, new_lab_matrix, estim_prob_labels, estim_prob_matrix, estim_LLH, computation_time, n_iter, status)     
    
end