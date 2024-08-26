using LinearAlgebra, Plots, Random, Distributions, DelimitedFiles, Statistics, PyCall, Clustering

include("Auxiliaries.jl")
include("SBM.jl")
include("LabeledGraph.jl")
include("SBMEstimator.jl")
include("Experiment.jl")

# --- Simulation parameters ----------------------------------------------------------------------------

n_simu = 1000
max_iter = 100
max_time = 10.0
form_prob_labels(x) = [(1-x)/3,1/3,(1+x)/3]
form_prob_matrix = "symmetric"
level = nothing
seed = nothing
simu_params = SimuParamSet(n_simu, form_prob_labels, form_prob_matrix, level, seed)

path_sims = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE2sym/Simulations"
path_estims = nothing
path_results = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE2sym"

# ------------------------------------------------------------------------------------------------------

# --- EXPE_identifiability : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(100)
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Sp+LloydL1", "Sp+LloydHuber0.05", "Sp+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Sp+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe2A = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_heterogeneity_identifiable : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(100)
param_prob_labels = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
p_intra = Parameter(0.9)
p_inter = Parameter(0.6)
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Sp+LloydL1", "Sp+LloydHuber0.05", "Sp+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Sp+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe2B1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_heterogeneity_non_identifiable : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(100)
param_prob_labels = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
p_intra = Parameter(0.9)
p_inter = Parameter(0.8)
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Sp+LloydL1", "Sp+LloydHuber0.05", "Sp+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Sp+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe2B2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_table : COMPARAISON DES METHODE Spectral+LLOYDHUBER TABLEAU RECAPITULATIF --------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(100)
param_prob_labels = Parameter([0.1,0.4,0.7])
p_intra = Parameter(0.9)
p_inter = Parameter([0.1, 0.6, 0.8])
noise_ini = Parameter(1.0)

output = "Table"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Sp+LloydL1", "Sp+LloydHuber0.05", "Sp+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Sp+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Sp+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe2C = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

expe_params = Dict{String,ExpeParamSet}("expe2A" => expe2A, "expe2B1" => expe2B1, "expe2B2" => expe2B2, "expe2C" => expe2C)
                                        
keys_expe = ["expe2A","expe2B1","expe2B2","expe2C"]

Experiment(simu_params, expe_params, keys_expe, path_results, path_sims, path_estims)