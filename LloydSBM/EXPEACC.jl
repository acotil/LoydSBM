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

path_sims = "/home/cotil/Documents/code_julia/SBM/Experiments/Work2/EXPE_ACC_sym/Simulations"
path_estims = nothing
path_results = "/home/cotil/Documents/code_julia/SBM/Experiments/Work2/EXPE_ACC_sym"

# ------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_table : COMPARAISON DES METHODE LLOYDHUBER TABLEAU RECAPITULATIF --------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter([10,30,50])
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.15 0.25 ; 0.35 0.45 ; 0.55 0.65 ; 0.75 0.85], ["0.2", "0.4", "0.6", "0.8"])
p_inter = Parameter([0.15 0.25 ; 0.35 0.45 ; 0.55 0.65 ; 0.75 0.85], ["0.2", "0.4", "0.6", "0.8"])
noise_ini = Parameter(1.0)

output = "Table"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["LloydL1","LloydMLE","GD"]
estimators["LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, with_ini = true)

benchmarks = ["accuracy","time"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe1_RI_sym_table = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_acc : RI + symmetric + accuracy ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["LloydL1","LloydMLE","GD"]
estimators["LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, with_ini = true)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe1_RI_sym_acc = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_n_lab : RI + symmetric + n_lab ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["LloydL1","LloydMLE","GD"]
estimators["LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, with_ini = true)

benchmarks = ["n_lab"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe1_RI_sym_n_lab = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_time : RI + symmetric + time ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["LloydL1","LloydMLE","GD"]
estimators["LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, with_ini = true)

benchmarks = ["time"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe1_RI_sym_time = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# -----------------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_iter : RI + symmetric + iter ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
noise_ini = Parameter(1.0)

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["LloydL1","LloydMLE","GD"]
estimators["LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, with_ini = true)
estimators["GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, with_ini = true)

benchmarks = ["iter"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K", "N", "h", "a", "b", "ω"]
expe1_RI_sym_iter = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# -----------------------------------------------------------------------------------------------------------------

# --- expe1_SI_sym_table : COMPARAISON DES METHODE LLOYDHUBER TABLEAU RECAPITULATIF --------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter([10,30,50])
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.15 0.25 ; 0.35 0.45 ; 0.55 0.65 ; 0.75 0.85], ["0.2", "0.4", "0.6", "0.8"])
p_inter = Parameter([0.15 0.25 ; 0.35 0.45 ; 0.55 0.65 ; 0.75 0.85], ["0.2", "0.4", "0.6", "0.8"])
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "Table"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1","Spectral+LloydMLE","Spectral+GD"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydMLE"] = SBMEstimator("LloydMLE", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+GD"] = SBMEstimator("GD", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy","time"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe1_SI_sym_table = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- expe1_SI_sym_acc : SI + symmetric + accuracy ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
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
expe1_SI_sym_acc = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_n_lab : RI + symmetric + n_lab ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
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

benchmarks = ["n_lab"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe1_SI_sym_n_lab = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_time : RI + symmetric + time ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
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

benchmarks = ["time"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe1_SI_sym_time = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# -----------------------------------------------------------------------------------------------------------------

# --- expe1_RI_sym_iter : RI + symmetric + iter ----------------------------------------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(30)
param_prob_labels = Parameter(0.0)
p_intra = Parameter([0.75 0.85],"0.8")
p_inter = Parameter([0.05 0.15 ; 0.15 0.25 ; 0.25 0.35 ; 0.35 0.45 ; 0.45 0.55 ; 0.55 0.65 ; 0.65 0.75 ; 0.75 0.85 ; 0.85 0.95], ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"])
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

benchmarks = ["iter"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K", "N", "h", "a", "b", "ω"]
expe1_SI_sym_iter = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# -----------------------------------------------------------------------------------------------------------------

expe_params = Dict{String,ExpeParamSet}("expe1_RI_sym_acc" => expe1_RI_sym_acc, 
                                        "expe1_RI_sym_time" => expe1_RI_sym_time,
                                        "expe1_RI_sym_n_lab" => expe1_RI_sym_n_lab,
                                        "expe1_RI_sym_iter" => expe1_RI_sym_iter,
                                        "expe1_RI_sym_table" => expe1_RI_sym_table,
                                        "expe1_SI_sym_acc" => expe1_SI_sym_acc,
                                        "expe1_SI_sym_n_lab" => expe1_SI_sym_n_lab,
                                        "expe1_SI_sym_time" => expe1_SI_sym_time,
                                        "expe1_SI_sym_iter" => expe1_SI_sym_iter,
                                        "expe1_SI_sym_table" => expe1_SI_sym_table)
                                         
keys_expe = ["expe1_RI_sym_acc", "expe1_RI_sym_n_lab", "expe1_RI_sym_time", "expe1_RI_sym_iter", "expe1_RI_sym_table", "expe1_SI_sym_acc", "expe1_SI_sym_n_lab", "expe1_SI_sym_time", "expe1_SI_sym_iter", "expe1_SI_sym_table"]

Experiment(simu_params, expe_params, keys_expe, path_results, path_sims, path_estims)
