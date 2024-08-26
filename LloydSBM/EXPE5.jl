using LinearAlgebra, Plots, Random, Distributions, DelimitedFiles, Statistics, PyCall, Clustering

include("Auxiliaries.jl")
include("SBM.jl")
include("LabeledGraph.jl")
include("SBMEstimator.jl")
include("Experiment.jl")

# --- Simulation parameters ----------------------------------------------------------------------------

n_simu = 1000
max_iter = 300
max_time = 10.0
form_prob_labels(x) = [(1-x)/3,1/3,(1+x)/3]
form_prob_matrix = "symmetric"
level = nothing
seed = nothing
simu_params = SimuParamSet(n_simu, form_prob_labels, form_prob_matrix, level, seed)

path_sims = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE5sym/Simulations"
path_estims = nothing
path_results = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE5sym"

# --- EXPE5A1 : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5A1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE1A1 : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["is_cv"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5A2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_heterogeneity_identifiable : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
p_intra = Parameter(0.9)
p_inter = Parameter(0.6)
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5B1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_heterogeneity_non_identifiable : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
p_intra = Parameter(0.9)
p_inter = Parameter(0.8)
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5B2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_table : COMPARAISON DES METHODE Spectral+LLOYDHUBER TABLEAU RECAPITULATIF --------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter([0.1,0.4,0.7])
p_intra = Parameter(0.9)
p_inter = Parameter([0.1, 0.6, 0.8])
noise_ini = Parameter(1.0)

output = "Table"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy","is_cv"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe5C = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

expe_params = Dict{String,ExpeParamSet}("expe5A1" => expe5A1, "expe5A2" => expe5A2, "expe5B1" => expe5B1, "expe5B2" => expe5B2, "expe5C" => expe5C)
                                         
keys_expe = ["expe5A1", "expe5A2", "expe5B1", "expe5B2", "expe5C"]

Experiment(simu_params, expe_params, keys_expe, path_results, path_sims, path_estims)


# --- Simulation parameters ----------------------------------------------------------------------------

n_simu = 1000
max_iter = 300
max_time = 10.0
form_prob_labels(x) = [(1-x)/3,1/3,(1+x)/3]
form_prob_matrix = "asymmetric"
level = nothing
seed = nothing
simu_params = SimuParamSet(n_simu, form_prob_labels, form_prob_matrix, level, seed)

path_sims = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE5asym/Simulations"
path_estims = nothing
path_results = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE5asym"

# --- EXPE5A1 : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5A1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE1A1 : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["is_cv"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5A2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_heterogeneity_identifiable : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
p_intra = Parameter(0.9)
p_inter = Parameter(0.6)
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5B1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_heterogeneity_non_identifiable : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
p_intra = Parameter(0.9)
p_inter = Parameter(0.8)
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe5B2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_table : COMPARAISON DES METHODE Spectral+LLOYDHUBER TABLEAU RECAPITULATIF --------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter([0.1,0.4,0.7])
p_intra = Parameter(0.9)
p_inter = Parameter([0.1, 0.6, 0.8])
noise_ini = Parameter(1.0)

output = "Table"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+VEM"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+VEM"] = SBMEstimator("VEM", tol=1e-7, tolFP=1e-5, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy","is_cv"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe5C = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

expe_params = Dict{String,ExpeParamSet}("expe5A1" => expe5A1, "expe5A2" => expe5A2, "expe5B1" => expe5B1, "expe5B2" => expe5B2, "expe5C" => expe5C)
                                         
keys_expe = ["expe5A1", "expe5A2", "expe5B1", "expe5B2", "expe5C"]

Experiment(simu_params, expe_params, keys_expe, path_results, path_sims, path_estims)
