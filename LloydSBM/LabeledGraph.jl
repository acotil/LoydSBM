"""
Structure representing the a graph whose each node have a label in [n_labels]

Arguments
----------
n_nodes : Int
    Total number of nodes

adj_matrix : Matrix{Int}
    The adjacency matrix of the graph. Each entry must belongs to {0,1}.

n_labels : Int
    Total number of labels.

lab_vector : Vector{Int}
    The Vector of labels. The i entry of this vector is the label of node i.

lab_matrix : Matrix{Float64}
    The matrix of labels, used for soft clustering. The (i,k) entry of this ùatrix correspond to the probability that the label of i being k.
"""
struct LabeledGraph
    n_nodes::Int                                # Number of nodes
    adj_matrix::Matrix{Int}                     # Adjacency matrix
    n_labels::Union{Nothing,Int}                # Number of labels
    lab_vector::Union{Nothing,Vector{Int}}          # Vector of labels
    lab_matrix::Union{Nothing,Matrix{Float64}}      # membership matrix of labels

    # Constructor function for struct LabGraph
    function LabeledGraph(adj_matrix::Matrix{Int}, lab_vector::Union{Nothing,Vector{Int}}, lab_matrix::Union{Nothing,Matrix{Float64}})::LabeledGraph

        # Check if adj_matrix is square matrix
        if size(adj_matrix,1) != size(adj_matrix,2)
            throw(ArgumentError("Matrix adj_matrix must be a square matrix."))
        end

        # Check adj_matrix 
        if ~all(in.(adj_matrix,([0,1],)))
            throw(ArgumentError("All elements of matrix adj_matrix must be 0 or 1."))
        end

        n_nodes = size(adj_matrix,1)

        if (lab_matrix !== nothing && lab_vector === nothing) || (lab_matrix === nothing && lab_vector !== nothing)
            throw(ArgumentError("If lab_vector !== nothing (resp. lab_matrix !== nothing) then lab_matrix !== nothing (resp. lab_vector !== nothing)"))
        elseif lab_vector !== nothing && lab_matrix !== nothing
            
            # Check if lab_vector and adj_matrix match
            if n_nodes!=size(lab_matrix,1) 
                throw(ArgumentError("sizes lab_matrix and adj_matrix do not match"))
            end
            if n_nodes!=length(lab_vector)
                throw(ArgumentError("sizes lab_vector and adj_matrix do not match"))
            end
            
            n_labels = size(lab_matrix,2)
            
            if lab_matrix_to_vector(lab_matrix) != lab_vector
                throw(ArgumentError("values of lab_vector and and lab_matrix do not match"))
            end
        else 
            n_labels = nothing
        end
    
        return new(n_nodes, adj_matrix, n_labels, lab_vector, lab_matrix)
    end
end

"""
Constructor of a LabeledGraph from an adjacency matrix and a vector of labels

Parameters
----------
adj_matrix : Matrix{Int}
    The adjacency matrix

lab_vector : Vector{Int}
    The vector of labels

Returns
-------
A new LabeledGraph. 
"""
function LabeledGraph(adj_matrix::Matrix{Int}, lab_vector::Vector{Int})::LabeledGraph
    return(LabeledGraph(adj_matrix, lab_vector, lab_vector_to_matrix(maximum(lab_vector), lab_vector)))
end

"""
Constructor of a LabeledGraph from an adjacency matrix and a matrix of labels

Parameters
----------
adj_matrix : Matrix{Int}
    The adjacency matrix

lab_matrix : Vector{Int}
    The matrix of labels

Returns
-------
A new LabeledGraph. 
"""
function LabeledGraph(adj_matrix::Matrix{Int}, lab_matrix::Matrix{Float64})::LabeledGraph
    return(LabeledGraph(adj_matrix, lab_matrix_to_vector(lab_matrix), lab_matrix))
end

"""
Constructor of a LabeledGraph from an adjacency matrix, set all the labels argument to nothing

Parameters
----------
adj_matrix : Matrix{Int}
    The adjacency matrix

Returns
-------
A new LabeledGraph. 
"""
function LabeledGraph(adj_matrix::Matrix{Int})
    return(LabeledGraph(adj_matrix, nothing, nothing))
end

"""
Constructor of a LabeledGraph. Generates a LabeledGraph from a Bernouilli SBM with self-loops.

Parameters
----------
n_nodes : Int
    The total number of nodes.

sbm : SBM
    A SBM structure.

Optional Arguments
----------
all_present : Bool
    If true then force the random initial labelisation to contain all the different labels

seed : Union{Nothing,Int}
    Seed for the random number generator.

Returns
----------
A new LabelGraph. 
"""
function LabeledGraph(n_nodes::Int, sbm::SBM ; all_present::Bool = false, seed::Union{Nothing,Int} = nothing)::LabeledGraph

    # Set the random seed if one was provided
    if seed !== nothing
        Random.seed!(seed)
    end

    # Construct the vector label
    if ~all_present
        lab_vector = [simprob(sbm.prob_labels) for i = 1:n_nodes]
    else
        while true
            lab_vector = [simprob(sbm.prob_labels) for i = 1:n_nodes]
            if length(Set(lab_vector))!=sbm.n_labels
                lab_vector = [simprob(sbm.prob_labels) for i = 1:n_nodes]
            else
                break
            end
        end
    end

    # Initialize adjacency matrix
    adj_matrix = [Int(rand() < sbm.prob_matrix[lab_vector[i],lab_vector[j]]) for i=1:n_nodes, j=1:n_nodes]

    n_labels,lab_vector = relabel(lab_vector)
    lab_matrix = lab_vector_to_matrix(n_labels,lab_vector)

    return(LabeledGraph(adj_matrix, lab_vector, lab_matrix))
end