#= COMMENTAIRES

#PRIORITAIRES

- Faire en sorte que la structure Parameter prenne en argument un type

- Recoder la structure ExpeParamSet pour mettre en remplaçant les argument par un nombre arbitraire de paramètres avec des noms différents

- modifier l'appelation des noms pour enlever le underscore dans latex

- Modifier randomize pour que le nouveau soit label soit différent du précédent

- Devra être ajouter l'argument "n_labels" à la structure Estimator car il faut différentier 
le nombre de label théorique (i.e le n_label de SBM), le nombre de label observer (i.e celui
de GraphLabeled) et le nombre de label désiré pour l'estimation (i.e celui de Estimator). 

- A chaque fois, ne pas oublier de modifier les description des fonctions.

- Ajouter initialisation clustering spectral et methode de merge en split 
(s'inspirer de sbm.R).

- Coder verbose et créer une structure results qui sera retournée par estimate si verbose = True

- FAIRE COMPARAISONS A sbm.R (à définir ici aussi car pas les mêmes que 
avant vu qu'il y a initialisation et merge and split en plus).

- Ajouter les types de sortie sur les fonctions (si possible).

- Ajouter dans les messages d'erreur relatifs à la validité d'une methode la plage des methode valide (exemple avec le status de IntParameter et FloatParameter) 

- mettre les arguments de toutes les fonctions en colonne (mise en page du code)

- Faire un abstract type Estimator et surcharger la fonction estimate

- Faire une super structure Experiments contenant simu_params, expe_params, keys_expe, path_results, path_sims, path_estims

#OPTIONNELS

- Dans Experiment, essayer d'avoir une meilleurs manière d'encoder la forme voulue de prob_labels (pour l'intant c'est une fonction form_prob_labels)

- Faire en sorte que estimate retourne une structure SBM

- Doit être modifiée la partie plot à la fin de main.jl pour qu'on ait pas à rentrer 
tout les noms des methodes à chaque fois.

- Dans les fonction Optimizers, mettre metric à la fin afin de mettre à la fin tout les 
arguments qui peuvent être égal à nothing.

- Pour aller plus vite, essayer de faire en sorte que la metric soit chargée avant utilisation
au même tire que la methode dans Estimate

=#


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
form_prob_matrix = "asymmetric"
level = nothing
seed = nothing
simu_params = SimuParamSet(n_simu, form_prob_labels, form_prob_matrix, level, seed)

path_sims = "/home/cotil/Documents/code_julia/SBM/Experiments/EXPE2asym/Simulations"
path_estims = nothing
path_results = "/home/cotil/Documents/code_julia/SBM/Experiments/EXPE2asym"

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
keys_estimators = ["Spectral","Spectral+LloydL1", "Spectral+LloydHuber0.025", "Spectral+LloydHuber0.05", "Spectral+LloydHuber0.075", "Spectral+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.025"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.025, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.075"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.075, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe_identifiability = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_heterogeneity_identifiable : COMPARAISON DES METHODES LLOYDHUBER PAR RAPPORT A L'IDENTIFIABILITE----------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(100)
param_prob_labels = Parameter([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
p_intra = Parameter(0.9)
p_inter = Parameter(0.7)
noise_ini = Parameter(1.0)
reg_coeff = 0.25

output = "GraphPB"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1", "Spectral+LloydHuber0.025", "Spectral+LloydHuber0.05", "Spectral+LloydHuber0.075", "Spectral+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.025"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.025, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.075"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.075, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe_heterogeneity_identifiable = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

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
keys_estimators = ["Spectral","Spectral+LloydL1", "Spectral+LloydHuber0.025", "Spectral+LloydHuber0.05", "Spectral+LloydHuber0.075", "Spectral+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.025"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.025, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.075"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.075, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]
expe_heterogeneity_non_identifiable = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

# --- EXPE_table : COMPARAISON DES METHODE Spectral+LLOYDHUBER TABLEAU RECAPITULATIF --------------------------------

# simulation parameters
n_labels = Parameter(3)
n_nodes = Parameter(100)
param_prob_labels = Parameter([0.1,0.4,0.7])
p_intra = Parameter(0.9)
p_inter = Parameter([0.6, 0.7, 0.8])
noise_ini = Parameter(1.0)

output = "Table"

# estimators parameters

estimators=Dict{String,SBMEstimator}()
keys_estimators = ["Spectral","Spectral+LloydL1", "Spectral+LloydHuber0.025", "Spectral+LloydHuber0.05", "Spectral+LloydHuber0.075", "Spectral+LloydL2"]
estimators["Spectral"] = SBMEstimator("Spectral", reg_coeff = reg_coeff)
estimators["Spectral+LloydL1"] = SBMEstimator("LloydMix", metric = "l1", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.025"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.025, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.05"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.05, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydHuber0.075"] = SBMEstimator("LloydMix", metric = "Huber", Huber_const = 0.075, max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)
estimators["Spectral+LloydL2"] = SBMEstimator("LloydMix", metric = "l2", max_iter = max_iter, max_time = max_time, reg_coeff = reg_coeff, with_ini = false)

benchmarks = ["accuracy"]

# building experiment
parameters = Dict{String,Parameter}(["K" => n_labels, 
                                    "N" => n_nodes, 
                                    "h" => param_prob_labels, 
                                    "a" => p_intra, 
                                    "b" => p_inter, 
                                    "ω" => noise_ini])
keys_parameters = ["K","N", "h", "a", "b", "ω"]

expe_table = ExpeParamSet(parameters, keys_parameters, estimators, keys_estimators, benchmarks, output)

# ------------------------------------------------------------------------------------------------------

expe_params = Dict{String,ExpeParamSet}("expe2Aasym" => expe_identifiability, "expe2B1asym" => expe_heterogeneity_identifiable, "expe2B2asym" => expe_heterogeneity_non_identifiable, "expe2Casym" => expe_table)
                                        
keys_expe = ["expe2Aasym","expe2B1asym","expe2B2asym","expe2Casym"]

Experiment(simu_params, expe_params, keys_expe, path_results, path_sims, path_estims)