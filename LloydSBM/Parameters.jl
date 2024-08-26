Status = ["DET_FIXED", "RAND_FIXED", "DET_VARYING", "RAND_VARYING"]

struct Parameter{T<:Real}
    status::String
    values::Dict{String,Vector{T}}
    keys_values::Vector{String}

    # Constructor function for struct IntParameter
    function Parameter(status::String, values::Dict{String,Vector{T}}, keys_values::Vector{String}) where {T<:Real}

        #Check if key_values are the keys of values 
        if ~all(haskey.([values],keys_values))
            throw(ArgumentError("keys_values must be the vector of keys of dictionnary values"))
        end
        

        #Check if values is in the right format
        if any(@. length(get([values],keys_values,"n")) != 2)
            throw(ArgumentError("the size of a Parameter values must be of length 2"))
        end

        # Check is the status is valid
        if status == "DET_FIXED"
            if length(values) != 1 
                throw(ArgumentError("If status of a Parameter is DET_FIXED then the size of values must be 1"))
            elseif any([values[k][1] != values[k][2] for k in keys_values])
                throw(ArgumentError("If status of a Parameter is DET_FIXED then the element of values must be of the form [x, x] where x::Real"))
            end
        elseif status == "RAND_FIXED"
            if length(values) != 1 
                throw(ArgumentError("If status of a Parameter is RAND_FIXED then the size of values must be 1"))
            elseif any([values[k][1] > values[k][2] for k in keys_values])
                throw(ArgumentError("If status of an Parameter is RAND_FIXED then the element of values must be of the form [x, y] where x <= y"))
            end
        elseif status == "DET_VARYING"
            if any([values[k][1] != values[k][2] for k in keys_values])
                throw(ArgumentError("If status of a Parameter is DET_VARYING then elements of values must be of the form [x, x] where x::Real"))
            end
        elseif status == "RAND_VARYING"
            if any([values[k][1] > values[k][2] for k in keys_values])
                throw(ArgumentError("If status of an Parameter is RAND_VARYING then elements of values must be of the form [x, y] where x <= y"))
            end
        else
            throw(ArgumentError("Invalid status: $status."))
        end
        return(new{T}(status, values, keys_values))
    end
end

function Parameter(values::Union{T,Vector{T},Matrix{T}}, keys_values::Union{Nothing,String,Vector{String}} = nothing) where {T<:Real}
    if typeof(values) == T
        if keys_values === nothing
            keys_values = string(values)
            return(Parameter("DET_FIXED", Dict(keys_values => [values, values]), [keys_values]))
        elseif typeof(keys_values) == String
            return(Parameter("DET_FIXED", Dict(keys_values => [values, values]), [keys_values]))
        else
            throw(ArgumentError("If a single values is given then a unique keys_values must be given"))
        end
    elseif typeof(values) == Vector{T}
        if keys_values === nothing
            keys_values = string.(values)
            return(Parameter("DET_VARYING", Dict([k => [v, v] for (v,k) in zip(values,keys_values)]), keys_values))
        elseif typeof(keys_values) == Vector{String}
            return(Parameter("DET_VARYING", Dict([k => [v, v] for (v,k) in zip(values,keys_values)]), keys_values))
        else
            throw(ArgumentError("If a DET_VARYING values is given then a vector keys_values must be given"))
        end
    elseif typeof(values) == Matrix{T}
        if size(values,1)==1
            if keys_values === nothing
                keys_values = string(round(500*(values[1,1]+values[1,2]))/1000)
                return(Parameter("RAND_FIXED", Dict(keys_values => [values[1,1], values[1,2]]), [keys_values]))
            elseif typeof(keys_values) == String
                return(Parameter("RAND_FIXED", Dict(keys_values => [values[1,1], values[1,2]]), [keys_values]))
            else
                throw(ArgumentError("If a single row matrix is given then a unique keys_values must be given"))
            end
        else
            if keys_values === nothing
                keys_values = string.(round(500*(values[:,1]+values[:,2]))/1000)
                return(Parameter("RAND_VARYING", Dict([keys_values[i] => [values[i,1], values[i,2]] for i in eachindex(keys_values)]), keys_values))
            elseif typeof(keys_values) == Vector{String}
                return(Parameter("RAND_VARYING", Dict([keys_values[i] => [values[i,1], values[i,2]] for i in eachindex(keys_values)]), keys_values))
            else
                throw(ArgumentError("If a RAND_VARYING values is given then a vector keys_values must be given"))
            end
        end
    end
end

function GetStatus(param::Parameter)
    return(param.status)
end

function GetName(param::Parameter)
    return(param.keys_values)
end

function GetValues(param::Parameter)
    if param.status == "DET_FIXED"
        return(param.values[param.keys_values[1]][1])
    elseif param.status == "RAND_FIXED"
        return(param.keys_values[1])
    elseif param.status == "DET_VARYING"
        return([param.values[k][1] for k in param.keys_values])
    elseif param.status == "RAND_VARYING"
        return(reduce(vcat,transpose.([param.values[k] for k in param.keys_values])))
    else
        throw(ArgumentError("invalid status in GetValue : $(param.status)"))
    end
end


