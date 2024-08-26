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

path_sims = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE4sym/Simulations"
path_estims = nothing
path_results = "/home/cotil/Documents/code_julia/SBM/Experiments/Work1/EXPE4sym"

# ------------------------------------------------------------------------------------------------------

# --- EXPE1A1 : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe4A1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE1A1 : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(50)
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["n_lab"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe4A2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

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
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe4B1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

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
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe4B2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE1D1 : COMPARAISON DES METHODE LLOYDHUBER PAR RAPPORT A LA TAILLE DU GRAPH (metric = accuracy) -----------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter([10,30,50])
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter(0.1)
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe4C1 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE1D2 : COMPARAISON DES METHODE LLOYDHUBER PAR RAPPORT A LA TAILLE DU GRAPH (metric = time) -----------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter([10,30,50])
param_prob_labels = Parameter(0.0)
p_intra = Parameter(0.9)
p_inter = Parameter(0.1)
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["time"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe4C2 = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

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
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe4D = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

expe_params = Dict{String,ExpeParamSet}("expe4A1" => expe4A1, "expe4A2" => expe4A2, "expe4B1" => expe4B1, "expe4B2" => expe4B2, "expe4C1" => expe4C1, "expe4C2" => expe4C2, "expe4D" => expe4D)
                                        
keys_expe = ["expe4A1","expe4A2","expe4B1","expe4B2","expe4C1","expe4C2","expe4D"]

Experiment(simu_params, expe_params, keys_expe, path_results, path_sims, path_estims)

