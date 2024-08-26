"""
Structure representing the law of a Bernouilli SBM with self-loops

Arguments
----------
n_labels : Int
    Total number of labels.

prob_labels : Vector{Float64}

    Probability vector representing the probability that a node has a certain label.
prob_matrix : Matrix{Float64}

    A matrix whose the (k,l) entry reperesents the probability to observe an edge from a node of label k to a node of label l.
"""
struct SBM
    n_labels::Int                       # Number of labels
    prob_labels::Vector{Float64}        # Probabilities that a node belongs to a community.
    prob_matrix::Matrix{Float64}        # Matrix of probabilities 

    # Constructor function for struct SBM
    function SBM(prob_labels::Vector{Float64}, prob_matrix::Matrix{Float64})

        # Check if prob_labels is a probability vector
        if ~(all(prob_labels .> 0) && 1.0-1/1000 < sum(prob_labels) < 1.0+1/1000)  
            throw(ArgumentError(string("Invalid value in argument prob_labels: The coefficients must be strictly positive and their sum must equal 1.")))
        end

        # Check if prob_matrix is a square matrix
        if size(prob_matrix,1) != size(prob_matrix,2)
            throw(ArgumentError("Matrix prob_matrix must be a square matrix."))
        end

        # Check if elements in prob_matrix are positive and smaller than 1
        if ~all(0 .<= prob_matrix .<= 1)
            throw(ArgumentError("All elements in prob_matrix must be strictly positive and smaller than 1."))
        end

        n_labels = size(prob_matrix,1)

        # Check if prob_labels matches the dimensions of B
        if length(prob_labels) != n_labels
            throw(ArgumentError("Number of elements in prob_labels must match the number of labels in SBM"))
        end
        
        return new(n_labels, prob_labels, prob_matrix)
    end
end
