
function SpectralOpt(n_nodes::Int, n_labels::Int, adj_matrix::Matrix{Int}, reg_coeff::Union{Nothing,Float64} = nothing)
    
    start = time()
    
    # Regularisation of the adjacency matrix
    if reg_coeff === nothing
        reg_adj_matrix = transpose(adj_matrix)*adj_matrix
    else
        reg_adj_matrix = adj_matrix .+ reg_coeff*mean(adj_matrix)
        reg_adj_matrix = transpose(reg_adj_matrix)*reg_adj_matrix
    end

    # Laplacian of the regularized adjacency matrix
    deg = vec(sum(reg_adj_matrix,dims=2))
    lap_adj_matrix = zeros(Float64, n_nodes, n_nodes)
    for i in 1:n_nodes, j in 1:n_nodes
        if deg[i] != 0 && deg[j] != 0
            lap_adj_matrix[i,j] = reg_adj_matrix[i,j]/sqrt(deg[i]*deg[j])
        end
    end

    # Eigen decomposition of the Laplacian matrix and the matrix of eigenvectors associated with the k largest eigenvalues
    eig = eigen(lap_adj_matrix)
    eigvalues = abs.(eig.values)
    eigvectors = real.(eig.vectors)
    U = zeros(Float64, n_labels, n_nodes)
    for k in 1:n_labels
        j = argmax(eigvalues)
        U[k,:] = eigvectors[:,j]
        eigvalues[j] = -1
    end

    # Normalisation of each rows of Union
    for i in 1:n_nodes
        s = sqrt(sum(U[:,i].^2))
        if s > 0.0
            U[:,i]=U[:,i]./s
        end
    end

    #Kmeans on rows of U
    estim_lab_vector = assignments(kmeans(U, n_labels))
    estim_lab_matrix = lab_vector_to_matrix(n_labels, estim_lab_vector)
    estim_prob_labels = CalculateProbLabels(estim_lab_matrix)
    estim_prob_matrix = CalculateProbMatrix(estim_lab_matrix, adj_matrix)
    estim_LLH = CalculateLLH(estim_lab_vector, estim_prob_matrix, adj_matrix)

    return(estim_lab_vector, estim_lab_vector, estim_lab_matrix, estim_prob_labels, estim_prob_matrix, estim_LLH, time() - start, 0, "CONVERGENCE")
end

function SpectralPyOpt(n_nodes::Int, n_labels::Int, adj_matrix::Matrix{Int}, reg_coeff::Union{Nothing,Float64} = nothing)
    
    start = time()
    
    # Regularisation of the adjacency matrix
    if reg_coeff === nothing
        reg_adj_matrix = transpose(adj_matrix)*adj_matrix
    else
        reg_adj_matrix = adj_matrix .+ reg_coeff*mean(adj_matrix)
        reg_adj_matrix = transpose(reg_adj_matrix)*reg_adj_matrix
    end

    cluster = pyimport("sklearn.cluster")
    SpectralClustering = cluster.SpectralClustering(n_labels, n_init = 10, affinity = "precomputed", assign_labels = "cluster_qr")
    SpectralClustering.fit(reg_adj_matrix)
    
    estim_lab_vector = Vector{Int64}(SpectralClustering.labels_ .+ 1)
    estim_lab_matrix = lab_vector_to_matrix(n_labels, estim_lab_vector)
    estim_prob_labels = CalculateProbLabels(estim_lab_matrix)
    estim_prob_matrix = CalculateProbMatrix(estim_lab_matrix, adj_matrix)
    estim_LLH = CalculateLLH(estim_lab_vector, estim_prob_matrix, adj_matrix)

    return(estim_lab_vector, estim_lab_vector, estim_lab_matrix, estim_prob_labels, estim_prob_matrix, estim_LLH, time() - start, 0, "CONVERGENCE")
end