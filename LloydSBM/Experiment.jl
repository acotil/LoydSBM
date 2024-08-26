include("SimuParamSet.jl")
include("ExpeParamSet.jl")

function Experiment(simu_params::SimuParamSet,
                    expe_params::Dict{String,ExpeParamSet},
                    keys_expe::Vector{String},
                    path_results::Union{Nothing,String}, 
                    path_sims::Union{Nothing,String} = nothing, 
                    path_estims::Union{Nothing,String} = nothing)
    
    #######################
    ## INITIALISATION

    # --- Loading parameters -------------------------------------------------------------------

    #Loading simulation parameters
    n_simu = simu_params.n_simu
    form_prob_labels = simu_params.form_prob_labels
    form_prob_matrix = simu_params.form_prob_matrix
    confidence_level = simu_params.confidence_level
    seed = simu_params.seed

    # Set the random seed if one was provided
    if seed !== nothing
        Random.seed!(seed)
    else
        seed = abs(rand(Int))
        Random.seed!(seed)
    end

    # ------------------------------------------------------------------------------------------

    # --- Creating basic dictionaries ----------------------------------------------------------

    expe_to_id = Dict{String,Vector{String}}([kexp => [""] for kexp in keys_expe])
    id_to_values = Dict{String,Matrix{Float64}}()
    id_to_expe = Dict{String,Vector{String}}()

    for (kexp,expe) in expe_params
        parameters = expe.parameters
        keys_parameters = expe.keys_parameters
        for kp in keys_parameters
            new_id = Vector{String}([])
            for id in expe_to_id[kexp], kval in parameters[kp].keys_values
                append!(new_id,[string(id,"$(kval)_")])
            end
            expe_to_id[kexp] = new_id
        end
        expe_to_id[kexp] = chop.(expe_to_id[kexp])
        for id in expe_to_id[kexp]
            values = reduce(vcat,transpose.([parameters[p].values[v] for (p,v) in zip(keys_parameters,split(id,"_"))]))
            if ~(id in keys(id_to_values))
                id_to_values[id] = values
                id_to_expe[id] = [kexp]
            elseif id_to_values[id] != values
                throw(ArgumentError("In all the given experiments, a key of a value of a parameter must refer to a unique value"))
            else
                append!(id_to_expe[id],[kexp])
            end
        end
    end

    n_total = n_simu*length(id_to_values)

    expe_to_results = Dict{String,Dict{Tuple{String,String,String},Vector{Float64}}}()
    for (kexp,expe) in expe_params
        expe_to_results[kexp] = Dict{Tuple{String,String,String},Vector{Float64}}()
        for id in expe_to_id[kexp], kben in expe.benchmarks, kest in expe.keys_estimators
            expe_to_results[kexp][id,kben,kest] = []
        end
    end

    # ------------------------------------------------------------------------------------------

    # --- Creating basic objects ---------------------------------------------------------------

    # If path_sims is given then we will save the simulations
    if path_sims === nothing
        save_sims = false
    else
        save_sims = true
        rm.(readdir(path_sims, join = true), recursive = true)  # If save_sims then clean all files in path_sims
    end

    # If path_sims is given then we will save the simulations
    if path_estims === nothing
        save_estims = false
    else
        save_estims = true
        rm.(readdir(path_estims, join = true), recursive = true)  # If save_estims then clean all files in path_estims
        
        for kexp in keys_expe
            path_exp = string(path_estims, "/$kexp")
            mkdir(path_exp)
        end
    end
    
    # ------------------------------------------------------------------------------------------

    # --- Saving simulation parameters ---------------------------------------------------------

    # If save_sims then creat a file which contains all the parameter of the simulation
    s = "GENERAL SIMULATION PARAMETERS :\n\n"
    s = string(s, "\tseed\n\t$seed\n\n")
    s = string(s, "\tn_simu\n\t$n_simu\n\n")
    s = string(s, "\tn_total\n\t$n_total\n\n")
    s = string(s, "\tform_prob_matrix\n\t$form_prob_matrix\n\n\n")
    for kexp in keys_expe
        expe = expe_params[kexp]
        s = string(s,"SIMULATIONS : $kexp\n\n")
        s = string(s, "\tn_param\n\t$(length(expe_to_id[kexp]))\n\n")
        for kpar in expe.keys_parameters
            param = expe.parameters[kpar]
            s = string(s, "\t$kpar\n\t$(GetStatus(param))\n\t$(GetName(param))\n\t$(GetValues(param))\n\n")
        end
        s = string(s, "\n")
    end

    if save_sims    
        open(string(path_sims, "/sims_params.txt"), "w") do f
            write(f, s)
        end
    else 
        open(string(path_results, "/sims_params.txt"), "w") do f
            write(f, s)
        end
    end
    
    # ------------------------------------------------------------------------------------------

    # --- Saving simulation parameters ---------------------------------------------------------
    s = ""
    for kexp in keys_expe
        expe = expe_params[kexp]
        for kesp in expe.keys_estimators
            estim = expe.estimators[kesp]
            s = string(s, "EXPERIMENT : $kexp ; ESTIMATOR \"$kesp\"\n\n")
            s = string(s, "\tmethod\n\t$(estim.method)\n\n")
            if estim.method in ["LloydForward","LloydBackward", "LloydMix","LloydMLE", "LloydTest", "GD"]
                s = string(s, "\tmax_iter\n\t$(estim.max_iter)\n\n")
                s = string(s, "\tmax_time\n\t$(estim.max_time)\n\n")
                s = string(s, "\tmetric\n\t$(estim.metric)\n\n")
                if estim.metric == "Huber"
                    s = string(s, "\tHuber_const\n\t$(estim.Huber_const)\n\n")
                end
                s = string(s, "\twith_ini\n\t$(estim.with_ini)\n\n")
                if ~(estim.with_ini)
                    if estim.reg_coeff === nothing
                        s = string(s, "\treg_coeff\n\t0\n\n\n")
                    else
                        s = string(s, "\treg_coeff\n\t$(estim.reg_coeff)\n\n\n")
                    end
                else
                    s = string(s, "\n")
                end
            elseif estim.method in ["VEM"]
                s = string(s, "\tmax_iter\n\t$(estim.max_iter)\n\n")
                s = string(s, "\tmax_time\n\t$(estim.max_time)\n\n")
                s = string(s, "\ttol\n\t$(estim.tol)\n\n")
                s = string(s, "\ttolFP\n\t$(estim.tolFP)\n\n")
                s = string(s, "\twith_ini\n\t$(estim.with_ini)\n\n")
                if ~(estim.with_ini)
                    if estim.reg_coeff === nothing
                        s = string(s, "\treg_coeff\n\t0\n\n\n")
                    else
                        s = string(s, "\treg_coeff\n\t$(estim.reg_coeff)\n\n\n")
                    end
                else
                    s = string(s, "\n")
                end
            elseif estim.method in ["Spectral","SpectralPy"]
                if estim.reg_coeff === nothing
                    s = string(s, "\treg_coeff\n\t0\n\n\n")
                else
                    s = string(s, "\treg_coeff\n\t$(estim.reg_coeff)\n\n\n")
                end
            end
        end
    end

    if save_estims
        open(string(path_estims, "/estims_params.txt"), "w") do f
            write(f, s)
        end
    else
        open(string(path_results, "/estims_params.txt"), "w") do f
            write(f, s)
        end
    end
    
    # ------------------------------------------------------------------------------------------

    #######################
    ## HEREDITY
    
    n_trial = 0

    for (id,values) in id_to_values, sim in 1:n_simu
        
        id_sim = string(id,":$sim")

        # --- Generating simulation parmeters --------------------------------------------------
        n_trial += 1
        println("$id_sim ; $n_trial/$n_total")
        
        n_labels = Int(floor(values[1,1] + rand() * (values[1,2] - values[1,1] + 1)))
        n_nodes = Int(floor(values[2,1] + rand() * (values[2,2] - values[2,1] + 1))) 
        
        param_prob_labels = values[3,1] + rand() * (values[3,2] - values[3,1])
        prob_labels = form_prob_labels(param_prob_labels)
        
        prob_matrix = zeros(n_labels,n_labels)
        
        p_intra = values[4,1] + rand() * (values[4,2] - values[4,1])
        p_inter = values[5,1] + rand() * (values[5,2] - values[5,1])
        if form_prob_matrix == "symmetric"
            for k in 1:n_labels, l in 1:n_labels
                if k==l
                   prob_matrix[k,l] = p_intra
                else
                    prob_matrix[k,l] = p_inter
                end
            end
        elseif form_prob_matrix == "asymmetric"
            for k in 1:n_labels, l in 1:n_labels
                if l == k
                   prob_matrix[k,l] = p_intra
                elseif l == (k+1)%n_labels
                    prob_matrix[k,l] = p_inter
                else
                    if p_inter < p_intra
                        prob_matrix[k,l] = p_inter+(p_intra-p_inter)/p_intra
                    else
                        prob_matrix[k,l] = p_inter+(p_intra-p_inter)/(1.0-p_intra)
                    end
                end
            end
        end
        
        noise_ini = values[6,1] + rand() * (values[6,2] - values[6,1])

        # --------------------------------------------------------------------------------------
        
        # --- Simulation -----------------------------------------------------------------------
    
        sbm = SBM(prob_labels, prob_matrix)

        true_graph = LabeledGraph(n_nodes, sbm, all_present = true)

        true_lab_vector = true_graph.lab_vector
        true_lab_matrix = true_graph.lab_matrix
        adj_matrix = true_graph.adj_matrix
        ini_lab_vector = randomized(n_labels, true_lab_vector, noise_ini, all_present = true)
        ini_lab_matrix = lab_vector_to_matrix(n_labels,ini_lab_vector)

        # --------------------------------------------------------------------------------------

        # --- Saving simulation ----------------------------------------------------------------

        if save_sims
            path_sim = string(path_sims, "/$id_sim")
            mkdir(path_sim)
            writedlm(string(path_sim, "/sim_params.txt"), ["n_labels" n_labels ; 
                                                        "n_nodes" n_nodes ; 
                                                        "param_prob_labels" param_prob_labels ; 
                                                        "p_intra" p_intra ; 
                                                        "p_inter" p_inter ; 
                                                        "noise_ini" noise_ini], '\t')
            writedlm(string(path_sim, "/prob_labels.txt"), prob_labels, '\t')
            writedlm(string(path_sim, "/prob_matrix.txt"), prob_matrix, '\t')
            writedlm(string(path_sim, "/true_lab_vector.txt"), true_lab_vector, '\t')
            writedlm(string(path_sim, "/adj_matrix.txt"), adj_matrix, '\t')
            writedlm(string(path_sim, "/ini_lab_vector.txt"), ini_lab_vector, '\t')
        end

        # --------------------------------------------------------------------------------------

        # --- ESTIMATION -------------------------------------------------------------------------
        
        for kexp in id_to_expe[id]

            for (kest,estim) in expe_params[kexp].estimators
        
                ini_lab_vector, estim_lab_vector, estim_lab_matrix, estim_prob_labels, estim_prob_matrix, estim_LLH, time, iter, estim_status = Estimate(estim, n_labels, adj_matrix, ini_lab_vector = ini_lab_vector, ini_lab_matrix = ini_lab_matrix)
                estim_accuracy = accuracy(true_lab_matrix, estim_lab_matrix)
                ini_accuracy = accuracy(true_lab_matrix, ini_lab_matrix)
                estim_RAG = RAG(ini_accuracy, estim_accuracy)
                estim_n_lab = Float64(length(Set(estim_lab_vector)))
                estim_results = Dict(["status" => estim_status, "iter" => iter, "is_cv" => float(estim_status == "CONVERGENCE"), "time" => time, "accuracy" => estim_accuracy, "RAG" => estim_RAG, "LLH" => estim_LLH, "n_lab" => estim_n_lab])
                for kben in expe_params[kexp].benchmarks
                    if (kest == "VEM" && kben != "is_cv")
                        if kben in keys(estim_results)
                            push!(expe_to_results[kexp][id,kben,kest],estim_results[kben]) 
                        else
                            throw(ArgumentError("The benchmark $kben given in experiment $kexp is not valid"))
                        end
                    else
                        if kben in keys(estim_results)
                            push!(expe_to_results[kexp][id,kben,kest],estim_results[kben])
                        else
                            throw(ArgumentError("The benchmark $kben given in experiment $kexp is not valid"))
                        end
                    end
                end

                if save_estims
                    path_estim = string(path_estims, "/$kexp/$id_sim:$kest")
                    mkdir(path_estim)
                    writedlm(string(path_estim, "/estim_results.txt"), estim_results, '\t')
                    writedlm(string(path_estim, "/estim_prob_labels.txt"), estim_prob_labels, '\t')
                    writedlm(string(path_estim, "/estim_prob_matrix.txt"), estim_prob_matrix, '\t')
                    writedlm(string(path_estim, "/estim_lab_vector.txt"), estim_lab_vector, '\t')
                    writedlm(string(path_estim, "/estim_lab_matrix.txt"), estim_lab_matrix, '\t')
                end

            end
        end
        
        # --------------------------------------------------------------------------------------
    end

    #######################
    ## RETURNS

    expe_to_mean = Dict{String,Dict{Tuple{String,String,String},Float64}}()
    
    if confidence_level !== nothing
        expe_to_CI = Dict{String,Dict{Tuple{String,String,String},Float64}}()
        norm = Normal()
        quant = quantile(norm, (1+confidence_level)/2)       
    else
        expe_to_CI = Dict{String,Nothing}()
    end
    
    for kexp in keys(expe_params)
        expe_to_mean[kexp] = Dict{Tuple{String,String,String},Float64}()
        if confidence_level !== nothing
            expe_to_CI[kexp] = Dict{Tuple{String,String,String},Float64}()
        else
            expe_to_CI[kexp] = nothing
        end
    end
    
    for (kexp,expe) in expe_params
        for k in keys(expe_to_results[kexp])
            expe_to_mean[kexp][k] = mean(expe_to_results[kexp][k])
            if confidence_level !== nothing
                expe_to_CI[kexp][k] = quant*stdm(expe_to_results[kexp][k], expe_to_mean[kexp][k])/sqrt(n_simu)
            end
        end
        Display(expe, kexp, expe_to_mean[kexp], expe_to_CI[kexp], path_results)
    end
    
end
