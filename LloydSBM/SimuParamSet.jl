struct SimuParamSet
    n_simu::Int
    form_prob_labels::Function
    form_prob_matrix::String
    confidence_level::Union{Nothing,Float64}
    seed::Union{Nothing,Int}
    
    function SimuParamSet(n_simu::Int,
                        form_prob_labels::Function,
                        form_prob_matrix::String,
                        confidence_level::Union{Nothing,Float64} = nothing,
                        seed::Union{Nothing,Int} = nothing)

        return(new(n_simu, form_prob_labels, form_prob_matrix, confidence_level, seed))
    end
end