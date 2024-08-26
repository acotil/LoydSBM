include("Parameters.jl")

Benchmarks = ["iter", "time", "accuracy", "LLH", "n_lab"]

struct ExpeParamSet
    parameters::Dict{String,Parameter}
    keys_parameters::Vector{String}
    estimators::Dict{String,SBMEstimator} 
    keys_estimators::Vector{String}
    benchmarks::Vector{String}
    output::Union{Nothing,String}
    #=
    Ajouter les vérifications suivante :
        1) vérifier n_labels positif
        2) vérifier n_nodes positif 
        3) vérifier que le type d'entré de form_pro_matrix est Float64 et que sont type de sortie est vecteur de proba de taille n_labels
        4) vérifier que la taille de form_prob_matrix est de taille n_labels x n_labels, que ses entrés sont dans {0,1} et que sa diagonal contient que des 1
        5) vérifier que p_intra et dans [0,1]
        6) idem pour p_inter
        7) idem pour noise_ini
        8) vérifier n_simu positif
        9) vérifier que les benchmarks appartiennent tous à Benchmarks
        10) si seed est donné, vérifier que c'est conforme (rechercher ce qui est attendu)
    =#
    function ExpeParamSet(parameters::Dict{String,Parameter},
                        keys_parameters::Vector{String},
                        estimators::Dict{String,SBMEstimator}, 
                        keys_estimators::Vector{String},
                        benchmarks::Vector{String},
                        output::Union{Nothing,String} = nothing)

        return(new(parameters,keys_parameters,estimators,keys_estimators,benchmarks,output))
    end
end

function DisplayTable(expe::ExpeParamSet, name_expe::String, expe_to_results::Dict{Tuple{String,String,String},Float64}, path_results::String)

    benchmarks = expe.benchmarks
    estimators = expe.keys_estimators
    varying_params = Vector{Parameter}([])
    keys_varying_params = Vector{String}([])
    varying_params_indicator = Vector{Bool}([])

    for kpar in expe.keys_parameters
        if expe.parameters[kpar].status in ["DET_VARYING","RAND_VARYING"]
            push!(varying_params,expe.parameters[kpar])
            push!(keys_varying_params,kpar)
            push!(varying_params_indicator,true)
        else
            push!(varying_params_indicator,false)
        end
    end

    results = Dict{Tuple{Vector{String},String,String},Float64}()
    for ((id,kben,kest),val) in expe_to_results
        results[split(id,"_")[varying_params_indicator],kben,kest] = round(val*1000)/1000
    end

    values = [[]]
    for param in varying_params
        copy_values = copy(values)
        values = []
        for v in copy_values
            for p in param.keys_values
                push!(values,vcat(v,p))
            end
        end
    end
    values = Vector{Vector{String}}(values)

    n_varying_params = size(varying_params)[1]
    n_benchmarks = size(benchmarks)[1]
    n_estimators = size(estimators)[1]
    
    lastline = "\\end{array}\n\$\$"
    firstline = chop("\$\$\n\\begin{array}{"*repeat("c",n_varying_params)*"|"*repeat(repeat("c",n_estimators)*"|",n_benchmarks), tail = 1)*"} \n"
    nameline1 = "\t"*chop(string((string.(keys_varying_params," & "))...) * string((string.(benchmarks," "*repeat("& ", n_estimators)))...),tail = 2)*"\\\\\n"
    nameline2 = "\t "*repeat("& ",n_varying_params)*chop(repeat(string(string.(estimators," & ")...),n_benchmarks), tail = 2)*"\\\\\n"
    hline = "\\hline\n"

    mainlines = ""
    previous_vals = ["" for v in 1:n_varying_params]
 
    for vals in values
        mainlines = string(mainlines,"\t")
        for v in 1:n_varying_params
            if vals[v] != previous_vals[v]
                if v == 1
                    mainlines = string(mainlines,hline)
                    mainlines = string(mainlines,"\t")
                end
                mainlines = string(mainlines,"$(vals[v]) & ")
            else 
                mainlines = string(mainlines," & ")
            end
        end
        for bench in benchmarks, estim in estimators
            mainlines = string(mainlines,"$(results[vals,bench,estim]) & ")
        end
        mainlines = chop(mainlines,tail = 2)
        mainlines = string(mainlines,"\\\\\n")
        previous_vals = copy(vals)
    end

    array = string(firstline,nameline1,nameline2,mainlines,lastline)
    open(string(path_results, "/$name_expe.txt"), "w") do f
        write(f, array)
    end
    return(nothing)
end

function DisplayGraphPB(expe::ExpeParamSet, name_expe::String, means::Dict{Tuple{String,String,String},Float64}, cis::Union{Nothing,Dict{Tuple{String,String,String},Float64}}, path_results::String)
    
    Y_label = expe.benchmarks[1]
    Curve_labels = expe.keys_estimators
    varying_params = Vector{Parameter}([])
    keys_varying_params = Vector{String}([])
    varying_params_indicator = Vector{Bool}([])

    for kpar in expe.keys_parameters
        if expe.parameters[kpar].status in ["DET_VARYING","RAND_VARYING"]
            push!(varying_params,expe.parameters[kpar])
            push!(keys_varying_params,kpar)
            push!(varying_params_indicator,true)
        else
            push!(varying_params_indicator,false)
        end
    end
    
    results_means = Dict{Tuple{String,String},Float64}()
    for ((id,kben,kest),val) in means
        results_means[kest,string(split(id,"_")[varying_params_indicator][1])] = val
    end
    
    param = varying_params[1]
    X_label = keys_varying_params[1]
    X_keys = param.keys_values
    
    X_values = Vector{Float64}([])
    if param.status == "DET_VARYING"
        X_values = GetValues(param)
    else
        val = GetValues(param)
        X_values = [(val[i,1]+val[i,2])/2 for i in 1:length(X_keys)]
    end
    id = []
    Y_values = reduce(vcat,transpose.([[results_means[kest, id] for kest in Curve_labels] for id in X_keys]))

    if cis !== nothing
        results_cis = Dict{Tuple{String,String},Float64}()
        for ((id,kben,kest),val) in cis
            results_cis[kest,string(split(id,"_")[varying_params_indicator][1])] = val
        end
        Y_errors = reduce(vcat,transpose.([[results_cis[kest, id] for kest in Curve_labels] for id in X_keys]))
        plot(X_values, Y_values,
            yerror = Y_errors, 
            label = reshape(Curve_labels,1,length(Curve_labels)),
            xlabel = X_label, ylabel = Y_label,
            xticks = (X_values,X_keys))
    else
        plot(X_values, Y_values,
            label = reshape(Curve_labels,1,length(Curve_labels)),
            xlabel = X_label, ylabel = Y_label,
            xticks = (X_values,X_keys))
    end
    
    savefig(string(path_results,"/$name_expe.png"))
    
    return(nothing)
end

function Display(expe::ExpeParamSet, name_expe::String, means::Dict{Tuple{String,String,String},Float64}, cis::Union{Nothing,Dict{Tuple{String,String,String},Float64}}, path_results::String)
    if expe.output == "Table"
        DisplayTable(expe, name_expe, means, path_results)
    elseif expe.output == "GraphPB"
        DisplayGraphPB(expe, name_expe, means, cis, path_results)
    else
        throw(ArgumentError("the output \"$(expe.output)\" of experiment $name_expe is not valid"))
    end
end