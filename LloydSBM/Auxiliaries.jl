
"""
Simulates a discrete distribution given by a probability vector

Parameters
----------
alpha : Vector{Float64}
    Probability vector.

Returns
-------
A random Int belonging to 1:length(alpha) 
"""
function simprob(alpha::Vector{Float64})::Int

    # Check if alpha is a probability vector
    if ~(all(alpha .> 0) || sum(alpha) == 1)  
        throw(ArgumentError(string("Invalid value in argument alpha: The coefficients must be strictly positive and their sum must equal 1.")))
    end
    
    # Simulation by autocorrelation
    r=rand()
    i=1
    cum=alpha[1]
    while r>cum
        i+=1
        cum+=alpha[i]
    end
    return(i)
end

"""
Convert a vector of labels to a matrix of labels

Parameters
----------
n_labels : Int
    Total number of labels

lab_vector : Vector{Int}
    vector of labels

Returns
-------
The corresponding membership matrix
"""
function lab_vector_to_matrix(n_labels::Int, lab_vector::Vector{Int})::Matrix{Float64}
    return([Float64(l == k) for l in lab_vector, k in 1:n_labels])
end

"""
Convert a matrix of labels to a vector of labels. Also work for a soft clustering matrix

Parameters
----------
lab_matrix : Matrix{Float64}
    matrix of labels

Returns
-------
The corresponding vector of labels
"""
function lab_matrix_to_vector(lab_matrix::Matrix{Float64})::Vector{Int}
    return(argmax.(eachrow(lab_matrix)))
end

"""
Relabels a label vector so that its labels range from 1 to n_labels

Parameters
----------
lab_vector : Vector{Int}
    vector of labels

Returns
-------
new_n_labels : Int
    the new number of labels

relab_vector : Vector{Int}
    the relabelized vector and 
"""
function relabel(lab_vector::Vector{Int})::Tuple{Int,Vector{Int}}
    relab_vector = copy(lab_vector)
    values = Set(relab_vector)
    max = maximum(relab_vector)
    for v in 1:max
        v = max-v+1
        if ~(v in values)
            relab_vector = [i-Int(i>v) for i in relab_vector]
        end
    end
    return(length(values),relab_vector)
end



"""
Check if the labels of a vector of labels range form 1 to n_labels

Parameters
----------
lab_vector : Vector{Int}
    A vector of labels.

Returns
-------
A Bool which is true if all labels are present.
"""
function allpresent(lab_vector::Vector{Int})::Bool
    return(length(Set(lab_vector)) == maximum(lab_vector))
end

"""
Check if the labels of a vector of labels range form 1 to n_labels for a given n_labels

Parameters
----------
lab_vector : Vector{Int}
    A vector of labels.

n_labels : Int
    The expected number of labels

Returns
-------
A Bool which is true if all labels are present.
"""
function allpresent(lab_vector::Vector{Int}, n_labels::Int)::Bool
    return((length(Set(lab_vector)) == n_labels) && (maximum(lab_vector) == n_labels))
end

"""
Check if the labels of a matrix of labels range form 1 to n_labels

Parameters
----------
lab_matrix : Matrix{Float64}
    A matrix of labels.

Returns
-------
A Bool which is true if all labels are present.
"""
function allpresent(lab_matrix::Matrix{Float64})::Bool
    return(all(sum(lab_matrix,dims = 1).>0.))
end

"""
Compute the MLE of the probability vector defining the probability that a node has a certain label from a vector of labels

Parameters
----------
n_labels : Int
    total number of labels (may be diffrent of the number of labels of vector labels)

lab_vector : Vector{Int}
    a vector of labels

Returns
-------
estim_prob_labels : Vector{Float64}
    the probability vector defining the probability that a node has a certain label
"""
function CalculateProbLabels(n_labels::Int,lab_vector::Vector{Int})::Vector{Float64}
    estim_prob_labels = zeros(n_labels)
    for l in lab_vector
        estim_prob_labels[l]+=1
    end
    estim_prob_labels/=length(lab_vector)
    return(estim_prob_labels)
end

"""
Compute the MLE of the probability vector defining the probability that a node has a certain label from a matrix of labels

Parameters
----------
lab_matrix : Matrix{Float64}
    a matrix of labels

Returns
-------
estim_prob_labels : Vector{Float64}
    the probability vector defining the probability that a node has a certain label
"""
function CalculateProbLabels(lab_matrix::Matrix{Float64})::Vector{Float64}
    return(vec(mean(lab_matrix,dims = 1)))
end

"""
Compute the MLE of the probability matrix of a Bernouilli SBM with self-loops from a vector of labels and an adjacency matrix

Parameters
----------
n_labels : Int
    Total number of labels (may be diffrent of the number of labels of vector labels).

lab_vector : Vector{Int}
    A vector of labels.

adj_matrix : Matrix{Int}
    The adjacency matrix of a Bernouilli SBM with self-loops.

Returns
-------
estim_prob_matrix : Matrix{Float64}
    A probability matrix of a Bernouilli SBM with self-loops.
"""
function CalculateProbMatrix(n_labels::Int, lab_vector::Vector{Int}, adj_matrix::Matrix{Int})::Matrix{Float64}
    n_nodes = length(lab_vector)
    estim_prob_matrix = zeros(Float64,n_labels,n_labels)
    for k in 1:n_labels, l in 1:n_labels
        estim_prob_matrix[k,l] = mean([adj_matrix[i,j] for i in 1:n_nodes, j in 1:n_nodes if (lab_vector[i]==k && lab_vector[j]==l)])
    end
    return(estim_prob_matrix)
end

"""
Compute the MLE of the probability matrix of a Bernouilli SBM with self-loops from a matrix of labels and an adjacency matrix

Parameters
----------
lab_matrix : Matrix{Float64}
    A matrix of labels.

adj_matrix : Matrix{Int}
    The adjacency matrix of a Bernouilli SBM with self-loops.

Returns
-------
estim_prob_matrix : Matrix{Float64}
    A probability matrix of a Bernouilli SBM with self-loops.
"""
function CalculateProbMatrix(lab_matrix::Matrix{Float64},adj_matrix::Matrix{Int})::Matrix{Float64}
    n_nodes,n_labels = size(lab_matrix)
    estim_prob_matrix = zeros(Float64,n_labels,n_labels)
    for k in 1:n_labels, l in 1:n_labels
        num = 0.0
        den = 0.0
        for i in 1:n_nodes, j in 1:n_nodes
            if lab_matrix[i,k] != 0.0 && lab_matrix[j,l] != 0.0
                if adj_matrix[i,j] == 1
                    num += lab_matrix[i,k] * lab_matrix[j,l]
                end
                den += lab_matrix[i,k] * lab_matrix[j,l]
            end
        end
        if den != 0
            estim_prob_matrix[k,l] = num / den
        end
    end
    return(estim_prob_matrix)
end

"""
Compute the Log-likelyhood of a Bernouilli SBM with self-loops from a vector of labels, a probability matrix and an adjacency matrix

Parameters
----------
lab_vector : Vector{Int}
    A vector of labels.

prob_matrix : Matrix{Float64}
    A probability matrix of a Bernouilli SBM with self-loops.

adj_matrix : Matrix{Int}
    The adjacency matrix of a Bernouilli SBM with self-loops.

Returns
-------
LLH : Float64
    the value of the computed Log-likelyhood of Bernouilli SBM with self-loops.
"""
function CalculateLLH(lab_vector::Vector{Int}, prob_matrix::Matrix{Float64}, adj_matrix::Matrix{Int})::Float64
    LLH=0.0
    n_nodes = length(lab_vector)
    for i in 1:n_nodes, j in 1:n_nodes
        if adj_matrix[i,j] == 1
            LLH += log(prob_matrix[lab_vector[i],lab_vector[j]])
        elseif adj_matrix[i,j] == 0
            LLH += log(1 - prob_matrix[lab_vector[i],lab_vector[j]])
        end
    end
    return(LLH)
end

"""
Compute the Log-likelyhood of a Bernouilli SBM with self-loops from a matrix of labels, a probability matrix and an adjacency matrix

Parameters
----------
lab_matrix : Matrix{Float64}
    A matrix of labels.

prob_matrix : Matrix{Float64}
    A probability matrix of a Bernouilli SBM with self-loops.

adj_matrix : Matrix{Int}
    The adjacency matrix of a Bernouilli SBM with self-loops.

Returns
-------
LLH : Float64
    the value of the computed Log-likelyhood of Bernouilli SBM with self-loops.
"""
function CalculateLLH(lab_matrix::Matrix{Float64},prob_matrix::Matrix{Float64},adj_matrix::Matrix{Int})::Float64
    return(CalculateLLH(lab_matrix_to_vector(lab_matrix),prob_matrix,adj_matrix))
end

"""
Compute the variationnal Log-likelyhood of a Bernouilli SBM with self-loops from a matrix of labels, a probability matrix and an adjacency matrix

Parameters
----------
lab_matrix : Matrix{Float64}
    A matrix of labels.

prob_matrix : Matrix{Float64}
    A probability matrix of a Bernouilli SBM with self-loops.

adj_matrix : Matrix{Int}
    The adjacency matrix of a Bernouilli SBM with self-loops.

Returns
-------
VARLLH : Float64
    the value of the computed Log-likelyhood of Bernouilli SBM with self-loops.
"""
function CalculateVARLLH(lab_matrix::Matrix{Float64}, prob_matrix::Matrix{Float64}, adj_matrix::Matrix{Int})::Float64
    VARLLH = 0.0
    n_nodes, n_labels = size(lab_matrix)
    for i in 1:n_nodes, k in 1:n_labels
        if lab_matrix[i,p] != 0.0
            a = ln(prob_matrix[k]/lab_matrix[i,p])
            for j in 1:n_nodes, l in 1:n_labels
                if lab_matrix[j,q] != 0.0
                    if adj_matrix[i,j] == 1
                        a += lab_matrix[j,q] * log(prob_matrix[k,l])
                    else
                        a += lab_matrix[j,q] * log(1.0 - prob_matrix[k,l])
                    end
                end
            end
            VARLLH += lab_matrix[i,p] * a
        end
    end
    return(VARLLH)
end

"""
Randomize a given vector of labels at rate p

Parameters
----------
n_labels : Int
    Total number of labels (may be diffrent of the number of labels of vector labels).

lab_vector : Vector{Int}
    A vector of labels.

p : Float64
    the probability to change a label of a node into a different label drawn uniformly at random into [n_label]

Optional Parameters
----------
all_present : Bool
    If true then repeat randomization until all labels are present at least once

Returns
-------
new_lab_vector : Vector{Int}
    A new random vector of labels.

WARNING
-------
If n_labels is large and all_present = true then this function can be very slow
"""
function randomized(n_labels::Int, lab_vector::Vector{Int}, p::Float64 ; all_present::Bool = false)::Vector{Int}
    
    n_nodes=length(lab_vector)
    
    if all_present && n_labels > length(lab_vector)
        throw(DomainError("the labels cannot all be present because the size of lab_vector is smaller than the number of labels requested"))
    end
    
    while true
        new_lab_vector=copy(lab_vector)
        for i in 1:n_nodes
            if rand()<p
                new_lab_vector[i]=rand(1:n_labels)
            end
        end
        if length(Set(new_lab_vector)) >= n_labels || ~all_present
            return(new_lab_vector)
        end
    end
end

"""
Transform a hard clustering matrix of labels into a soft clustering matrix of labels by changing all the zero values into epsilon.

Parameters
----------
lab_matrix : Matrix{Float64}
    A matrix of labels.

epsilon : Float64
    the value by which the zero values will be changed

Returns
-------
new_lab_matrix : Matrix{Float64}
    A new random vector of labels.
"""
function soften(lab_matrix::Matrix{Float64}, epsilon::Float64)::Matrix{Float64}
    n_nodes, n_labels = size(lab_matrix)
    new_lab_matrix = lab_matrix .+ (epsilon/(n_labels*(1.0-epsilon)-1.0))
    for i in 1:n_nodes
        new_lab_matrix[i,:]./=sum(new_lab_matrix[i,:])
    end
    return(new_lab_matrix)
end

"""
Set the label of a given node to a given label

Parameters
----------
lab_vector : Vector{Int}
    A vector of labels.

i : Int
    The indice of the node

k : Int
    The new label of the given node

Returns
-------
new_labels : Vector{Int}
    A new vector of labels.
"""
function set(lab_vector::Vector{Int},i::Int,k::Int)::Vector{Int}
    new_lab_vector = copy(lab_vector)
    new_lab_vector[i] = k
    return(new_lab_vector)
end

"""
Compute the pseudo-Huber loss function of parameter a
Detail : pseudo-Huber(x,alpha) = (a^2)*sqrt(1+(x/a)^2)-a^2

Parameters
----------
x : Float64
    the argument of the pseudo-huber loss function
a : Float64
    the parameter of the pseudo-huber loss function

Returns
-------
The value of the pseudo-Huber loss function of parameter a at x
"""
function Huber(x::Float64,a::Float64)::Float64
    if abs(x)<=a
        return((x^2)/2)
    else
        return(a*(abs(x)-a/2))
    end
end

"""
Compute the metric for LloydForward and LloydBackward method

Parameters
----------
metric : String
    the name of the used metric

estim_pi : Vector{Float64}
    the MLE of the vector of probabilities to oberve an edge between node i and a node of label l

estim_prob_matrix : Vector{Float64}
    the MLE of the vector of probabilities to oberve an edge between a node of label k and a node of label l

x : Vector{Int}
    The observed number of edges between node i and nodes of label l

N : Vector{Int}
    The number of node of a given label l

Huber_const : Float64
    The parameter of the pseudo-Huber loss function if metric is Huber

pen : Float64
    The value of the variationnal penalisation

Returns
-------
The value of the metric for the given arguments and parameters
"""
function Dist(metric::String, estim_pi::Vector{Float64}, estim_prob_matrix::Vector{Float64}, x::Vector{Int}, N::Vector{Int}, Huber_const::Union{Float64,Nothing}=nothing, pen::Union{Float64,Nothing}=nothing)::Float64
    if metric == "l1" # Distance l1
        return(sum(abs.(estim_pi .- estim_prob_matrix)))
    elseif metric == "l2" # Distance l2
        return(sum((estim_pi .- estim_prob_matrix).^2))
    elseif metric == "Huber"
        return(sum(Huber.(estim_pi .- estim_prob_matrix, Huber_const)))
    elseif metric == "Ent" # Distance MLE sans la correction par le nombre de noeuds
        return(sum(.- estim_pi .* log.(estim_prob_matrix .+ Float64.(estim_prob_matrix.<=0.0)) .- (1 .- estim_pi) .* log.(1 .- estim_prob_matrix .+ Float64.(estim_prob_matrix.>=1.0))))
    elseif metric == "Rent" # Distance MLE sans la correction par le nombre de noeuds
        return(sum(.- estim_pi .* log.((estim_prob_matrix .+ Float64.(estim_prob_matrix.<=0.0))./(estim_pi.+ Float64.(estim_pi.<=0.0))) .- (1 .- estim_pi) .* log.((1 .- estim_prob_matrix .+ Float64.(estim_prob_matrix .>= 1.0))./(1 .- estim_pi.+ Float64.(estim_pi.>=1.0)))))
    elseif metric == "linfty" # Distance linfty
        return(maximum(abs.(estim_pi .- estim_prob_matrix)))
    elseif metric == "l1c" # Distance l1 avec la correction par le nombre de noeuds
        return(sum(abs.(x .- N .* estim_prob_matrix)))
    elseif metric == "l2c" # Distance l2 avec la correction par le nombre de noeuds
        return(sum((x .- N .* estim_prob_matrix).^2))            
    elseif metric == "linftyc" # Distance linfty avec la correction par le nombre de noeuds
        return maximum(abs.(x .- N .* estim_prob_matrix))
    elseif metric == "Entc" # Distance MLE
        return(sum(.- x .* log.(estim_prob_matrix .+ Float64.(estim_prob_matrix.<=0.0)) .- (N .- x) .* log.(1 .- estim_prob_matrix .+ Float64.(estim_prob_matrix.>=1.0))))            
    elseif metric == "l1p" # Distance l1
        return(sum(abs.(estim_pi .- estim_prob_matrix))-0.5*log(pen))
    elseif metric == "l2p" # Distance l2
        return(sum((estim_pi .- estim_prob_matrix).^2 .-log.(N/sum(N)).^2))
    elseif metric == "Huberp"
        return(sum(Huber.(estim_pi .- estim_prob_matrix, Huber_const).-log.(N/sum(N))))
    elseif metric == "Entp" # Distance MLE sans la correction par le nombre de noeuds
        return(sum(.- estim_pi .* log.(estim_prob_matrix .+ Float64.(estim_prob_matrix.<=0.0)) .- (1 .- estim_pi) .* log.(1 .- estim_prob_matrix .+ Float64.(estim_prob_matrix.>=1.0)).-log.(N/sum(N))))
    elseif metric == "Rentp" # Distance MLE sans la correction par le nombre de noeuds
        return(sum(.- estim_pi .* log.((estim_prob_matrix .+ Float64.(estim_prob_matrix<=0.0))./(estim_pi.+ Float64.(estim_pi.<=0.0))) .- (1 .- estim_pi) .* log.((1 .- estim_prob_matrix .+ Float64.(estim_prob_matrix .>= 1.0))./(1 .- estim_pi.+ Float64.(estim_pi.>=1.0))).-log.(N/sum(N))))
    elseif metric == "linftyp" # Distance linfty
        return(maximum(abs.(estim_pi .- estim_prob_matrix).-log.(N/sum(N))))
    elseif metric == "l1cp" # Distance l1 avec la correction par le nombre de noeuds
        return(sum(abs.(x .- N .* estim_prob_matrix).-log.(N/sum(N))))
    elseif metric == "l2cp" # Distance l2 avec la correction par le nombre de noeuds
        return(sum((x .- N .* estim_prob_matrix).^2 .-log.(N/sum(N))).^2)            
    elseif metric == "linftycp" # Distance linfty avec la correction par le nombre de noeuds
        return maximum(abs.(x .- N .* estim_prob_matrix.-log.(N/sum(N))))
    elseif metric == "Entcp" # Distance MLE
        return(sum(.- x .* log.(estim_prob_matrix .+ Float64.(estim_prob_matrix.<=0.0)) .- (N .- x) .* log.(1 .- estim_prob_matrix .+ Float64.(estim_prob_matrix.>=1.0)))-log(pen))            
    else
        throw(ArgumentError("invalid metric ($metric) given to Dist"))
    end
end

"""
Compute the value of the objective function of the Lloyd type algorithm for a Bernouilli SBM with self-loops from a vector of labels, a probability matrix and an adjacency matrix

Parameters
----------
metric : String
    the name of the used metric

n_nodes : Int
    the number of nodes

n_labels : Int
    Total number of labels (may be diffrent of the number of labels of vector labels).

lab_vector : Vector{Int}
    A vector of labels.

estim_pi : Vector{Float64}
    the MLE of the vector of probabilities to oberve an edge between node i and a node of label l

estim_prob_matrix : Vector{Float64}
    the MLE of the vector of probabilities to oberve an edge between a node of label k and a node of label l

x : Vector{Int}
    The observed number of edges between node i and nodes of label l

N : Vector{Int}
    The number of node of a given label l

Huber_const : Float64
    The parameter of the pseudo-Huber loss function if metric is Huber

Returns
-------
obj : Float64
    the value of the objective function of the Lloyd type algorithm for a Bernouilli SBM with self-loops
"""
function CalculateObjective(metric::String, n_labels::Int, n_nodes::Int, lab_vector::Vector{Int}, estim_pi::Matrix{Float64}, estim_prob_matrix::Matrix{Float64}, x::Matrix{Int}, N::Vector{Int}, Huber_const::Union{Float64,Nothing}=nothing)::Float64
    obj = 0.0
    for k in 1:n_labels
        for i in 1:n_nodes
            if lab_vector[i] == k
                obj += Dist(metric, estim_pi[i, :], estim_prob_matrix[k, :], x[i, :], N, Huber_const, Float64(N[k]/n_labels))
            end
        end
    end
    return(obj)
end

"""
Swap to element in a permutation vector

Parameters
----------
perm : Vector{Int}
    a permutation vector

i : Int
    the first element to swap

j : Int
    the second element to swap
"""
function Swap!(perm::Vector{Int}, i::Int, j::Int)
    current_element = perm[j]
    perm[j]=perm[i]
    perm[i] = current_element
end

"""
Compute the permutation following a given permutation vector in the lexical order

Parameters
----------
n_labels : Int
    Total number of labels (may be diffrent of the number of labels of vector labels).

perm : Union{Nothing,Vector{Int}}
    a permutation vector or nothing

Returns
-------
next_perm : Union{Nothing,Vector{Int}}
    the next permutation in the lexical order or nothing if there is no next permutation
"""
function NextPerm(n_labels::Int, perm::Union{Nothing,Vector{Int}})::Union{Nothing,Vector{Int}}
    if perm !== nothing
        next_perm = copy(perm)
        i = n_labels-1
        while i>0 && perm[i]>=perm[i+1]
            i-=1
        end
        if i == 0
            return(nothing)
        else
            j = n_labels
            while next_perm[j] < next_perm[i]
                j-=1
            end
            Swap!(next_perm,i,j)
            i += 1
            j = n_labels
            while j > i
                Swap!(next_perm,i,j)
                j-=1
                i+=1
            end
            return(next_perm)
        end
    else
        return(nothing)
    end
end

"""
Compute the trace of the matrix perm(mat)[i,j] = mat[i,perm[j]]. T is a subtype of Real

Parameters
----------
mat : Matrix{T}
    the matrix.

perm : Vector{Int}
    the permutation

n : Int
    the size of the matrix

Returns
-------
trace_perm : T
    the value of the trace of perm(mat)
"""
function TracePerm(mat::Matrix{T}, perm::Vector{Int}, n::Int)::T where {T<:Real}
    return(return(sum([mat[i,perm[i]] for i in 1:n])))
end

"""
Compute the largest trace of the matrix perm(mat)[i,j] = mat[i,perm[j]] for all permutation vector. T is a subtype of Real

Parameters
----------
mat : Matrix{T}
    the matrix.

n : Int
    the size of the matrix

Returns
-------
best_trace : T
    the value of largest trace of perm(mat) for all permutation vector

best_perm : Vector{Int}
    the corresponding permutation
"""
function MaxTrace(mat::Matrix{T},n::Int)::Tuple{T,Vector{Int}} where {T<:Real}
    current_perm = [i for i in 1:n]
    best_perm = copy(current_perm)
    best_trace = -Inf
    while current_perm !== nothing
        current_trace = TracePerm(mat,current_perm,n)
        if current_trace > best_trace
            best_trace = current_trace
            best_perm = copy(current_perm)
        end
        current_perm = NextPerm(n,current_perm)
    end
    return(best_trace, best_perm)
end

"""
Compute the rate of common label between two vectors of labels up to a permutation on labels.

Parameters
----------
n_labels : Int
    Total number of labels (may be diffrent of the number of labels of vector labels).

lab_vector1 : Vector{Int}
    the first vector of labels

lab_vector2 : Vector{Int}
    the second vector of labels

Returns
-------
the rate of common label between the two vectors of labels
"""
function accuracy(n_labels::Int, lab_vector1::Vector{Int}, lab_vector2::Vector{Int})::Float64
    A = zeros(Int, n_labels, n_labels)
    n_nodes=length(lab_vector1)
    for i in 1:n_nodes
        A[lab_vector1[i], lab_vector2[i]] += 1
    end
    max_trace = Float64(MaxTrace(A,n_labels)[1])
    return(max_trace/n_nodes)
end

"""
Test if two vectors of labels are equal up to a permutation on labels.

Parameters
----------
n_labels : Int
    Total number of labels (may be diffrent of the number of labels of vector labels).

lab_vector1 : Vector{Int}
    the first vector of labels

lab_vector2 : Vector{Int}
    the second vector of labels

Returns
-------
a Bool which is true if the two vectors of labels are equal up to a permutation
"""
function islabequal(n_labels::Int, lab_vector1::Vector{Int}, lab_vector2::Vector{Int})::Bool
    current_perm = [i for i in 1:n_labels]
    perm_vector2 = copy(lab_vector2)
    equal = true
    while ~isequal(lab_vector1, perm_vector2)
        current_perm = NextPerm(n_labels,current_perm)
        if current_perm === nothing 
            equal = false
            break
        end
        perm_vector2 = [current_perm[i] for i in perm_vector2]        
    end
    return(equal)
end

"""
Compute the rate of common label between two matrix of labels up to a permutation on labels.

Parameters
----------
lab_matrix1 : Matrix{Float64}
    the first  matrix of labels

lab_matrix2 : Matrix{Float64}
    the second matrix of labels

Returns
-------
the rate of common label between the two matrix of labels
"""
function accuracy(lab_matrix1::Matrix{Float64}, lab_matrix2::Matrix{Float64})::Float64
    n_nodes, n_labels=size(lab_matrix1)
    A = zeros(Float64, n_labels, n_labels)
    for k in 1:n_labels, l in 1:n_labels, i in 1:n_nodes
        A[k,l] += minimum([lab_matrix1[i,k],lab_matrix2[i,l]])
    end
    max_trace = Float64(MaxTrace(A,n_labels)[1])
    return(max_trace/n_nodes)
end

"""
Test if two vectors of labels are equal up to a permutation on labels and a tolerance threshold.

Parameters
----------
lab_matrix1 : Matrix{Float64}
    the first  matrix of labels

lab_matrix2 : Matrix{Float64}
    the second matrix of labels

tol : Float64
    the tolerance threshold (1e-5 by default)

Returns
-------
a Bool which is true if the two vectors of labels are equal up to a permutation and a tolerance threshold.
"""
function islabequal(lab_matrix1::Matrix{Float64}, lab_matrix2::Matrix{Float64}, tol::Float64 = 1e-5)::Bool
    n_nodes = size(lab_matrix1,1)
    return(accuracy(lab_matrix1,lab_matrix2) >= 1.0-tol/n_nodes)
end

function RAG(ini_accuracy::Float64, estim_accuracy::Float64)::Float64
    if estim_accuracy > ini_accuracy
        return((estim_accuracy - ini_accuracy)/(1.0-ini_accuracy))
    else
        return((estim_accuracy - ini_accuracy)/ini_accuracy)
    end
end