# As of 2024-09-20 This is a Monte Carlo simulation showing the Hardy-Weinberg principle using a "AA", "Aa" and "aa" genotype
# future interation to include any single gene letters 
using Pkg
Pkg.activate(".")

include("GeneString.jl")
using .GeneString
using CairoMakie
using Chain
using StatsBase




using Base.Iterators
using Random
using Statistics
using Test




gene = "XX"
freq = 0.5
GeneString.total_genome_single_gene(gene)
GeneString.compute_genotype_HW_frequency(freq)
generation = ["XX"]
GeneString.compute_genotype_counts_frequency(generation)
GeneString.produce_1_offspring(parent1,parent2)
GeneString.isHomozygous(gene)
GeneString.isHeterozygous(gene)
GeneString.isHomozygousDominant(gene)
GeneString.isHomozygousRecessive(gene)
GeneString.detectZygosity(gene)
genotype = "QqAABbCc"
gen_iter = GeneString.split_genotype_into_pairs_iterator(genotype)
GeneString.isHomozygous.(gen_iter) 
GeneString.isHomozygousDominant.(gen_iter)
GeneString.isHomozygousRecessive.(gen_iter)
GeneString.isHeterozygous.(gen_iter) 
GeneString.detectZygosity.(gen_iter)
genotypes = ["AA", "Aa", "aa"]
p = 0.5
population_size = 1_000
ran_off = GeneString.RandomNumberOffspring(2,4)
fix_off = GeneString.FixedNumberOffspring(1)
weights = GeneString.compute_genotype_HW_frequency(p)
generation = GeneString.sample_initial_population(genotypes,weights,population_size)
GeneString.number_offspring(ran_off)
GeneString.number_offspring(fix_off)
GeneString.next_generation(fix_off,generation)
GeneString.next_generation(ran_off,generation)


# Create just AA: OX allele and aa: wild type
num_gens = 20
genotypes = ["AA","aa"]
p = 0.5
population_size = 1000
weights = [p,1-p]
ran_off = GeneString.RandomNumberOffspring(2,4)
fix_off = GeneString.FixedNumberOffspring(1)

uncapped_gen = []
population = GeneString.sample_initial_population(genotypes,weights,population_size)
sex = ascribe_sex_based_on_gene(population)
inds = indices_to_keep(population,sex)
push!(uncapped_gen,(population,sex))


capped_gen = []
push!(capped_gen,(population[inds],sex[inds]))
  

function test_pop1!(uncapped_gen,num_gens,ran_off,population,sex)
    for i in Base.OneTo(num_gens)
        population = GeneString.next_generation(ran_off ,population,sex)
        sex = ascribe_sex_based_on_gene(population)
        push!(uncapped_gen,(population,sex))
    end
    uncapped_gen
end



function test_pop2!(capped_gen,num_gens,ran_off,population,sex)
    for i in Base.OneTo(num_gens)
        population = GeneString.next_generation(ran_off ,population,sex)
        sex = ascribe_sex_based_on_gene(population)
        inds = indices_to_keep(population,sex)
        population,sex = population[inds],sex[inds]
        push!(capped_gen,(population,sex))
    end
    capped_gen
end

test_pop1!(uncapped_gen,num_gens,ran_off,population,sex)
test_pop2!(capped_gen,num_gens,ran_off,population,sex)



capped_gen_counts = count_unique.([c[1] for c in capped_gen])
uncapped_gen_counts = count_unique.([c[1] for c in uncapped_gen])



pop_gen = @chain capped_gen_counts begin
    filter.(p -> first(p) != "AA" , _)
    values.(_)
    sum.(_)
end


uncapped_pop_gen = @chain uncapped_gen_counts begin
    filter.(p -> first(p) != "AA" , _)
    values.(_)
    sum.(_)
end

end_diff = Int(round(100*(uncapped_pop_gen[end]-pop_gen[end])/uncapped_pop_gen[end],digits = 0))
f = Figure(size = (800,800),fontsize = 22)
ax = Axis(f[1,1],xlabel = "Generation",ylabel = "Population",aspect = 1,title = "Female removal if born with Oxitec gene versus wild type",subtitle = "N₀ = $(population_size), Ox₀ ratio = $(p), Diff at end ≈ $(end_diff)%")
lines!(ax,Base.OneTo(num_gens+1),pop_gen,label="Removed female\nwith Ox gene")
lines!(ax,Base.OneTo(num_gens+1),uncapped_pop_gen,label="Wild type",linestyle = :dash)
axislegend("Population")
f
save("female_removed_with_ox_init_rate-$(p)_N0-$(population_size)_diff-$(end_diff).png",f)

count_unique(sex)
count_unique(population)


count_unique(sex)
count_unique(population)


population = GeneString.next_generation(ran_off ,population)
sex = ascribe_sex_based_on_gene(population)
count_unique(sex)
count_unique(population)
population,sex = remove_female_with_ox_gene(population,sex)
count_unique(sex)
count_unique(population)

population = GeneString.next_generation(ran_off ,population)
sex = ascribe_sex_based_on_gene(population)
population,sex = remove_female_with_ox_gene(population,sex)
count_unique(sex)
count_unique(population)


p = 0.5
weights = GeneString.compute_genotype_HW_frequency(p)
population_size_const = 1_000_000
population_size_grow = 100
num_gens = 15
number_sims = 1
genotypes = ["XX", "Xx", "xx"]
ran_off = GeneString.RandomNumberOffspring(1,3)
fix_off = GeneString.FixedNumberOffspring(1)



# Constant size
generations_const_pop = Vector{String}[]
N₀_const = GeneString.sample_initial_population(genotypes,weights,population_size_const)
push!(generations_const_pop,N₀_const)
epoch_const = GeneString.evolve_generations(fix_off,num_gens,generations_const_pop)
stats_per_epoch_const = GeneString.evolve_stats(epoch_const)


# Grow size
generations_grow_pop = Vector{String}[]
N₀_grow = GeneString.sample_initial_population(genotypes,weights,population_size_grow)
push!(generations_grow_pop,N₀_grow)
epoch_grow = GeneString.evolve_generations(ran_off,num_gens,generations_grow_pop)
stats_per_epoch_grow = GeneString.evolve_stats(epoch_grow)



f = Figure(size = (800,800),fontsize = 22)
ax = Axis(f[1,1],xlabel = "Generation",ylabel = "Population",title = "Population growth",aspect = 1)
lines!(ax,1:num_gens+1,length.(epoch_const), label="Constants") 
lines!(ax,1:num_gens+1,length.(epoch_grow), label="Growing")
axislegend("Population")
f
save("population_count_constant_exponential.png",f)


pconst = (
    genotypes = ["XX", "Xx", "xx"],
    num_sims = 10,
    num_generations = 15,
    population_size = 10_000, 
    p = 0.5,
    offspring = GeneString.FixedNumberOffspring(1))
weights = GeneString.compute_genotype_HW_frequency(pconst[:p])
sims_const = GeneString.simulate_multiple_epochs(
    pconst[:offspring],
    pconst[:num_generations],
    pconst[:num_sims],
    pconst[:genotypes],
    weights,
    pconst[:population_size])
    


pgrow = (
    genotypes = ["XX", "Xx", "xx"],
    num_sims = 10,
    num_generations = 15,
    population_size = 50, 
    p = 0.5,
    offspring = GeneString.RandomNumberOffspring(1,3))
weights = GeneString.compute_genotype_HW_frequency(pgrow[:p])
sims_grow = GeneString.simulate_multiple_epochs(
    pgrow[:offspring],
    pgrow[:num_generations],
    pgrow[:num_sims],
    pgrow[:genotypes],
    weights,
    pgrow[:population_size])

    




# Execute a single evolution
f = GeneString.plot_single_epoch_sample_path(pgrow,sims_const[1]["generation"])
save("monte_carlo_simulation_verify_hardy_weinberg_single_sim_const_pop.png",f)

f = GeneString.plot_single_epoch_sample_path(pgrow,sims_grow[1]["generation"])
save("monte_carlo_simulation_verify_hardy_weinberg_single_sim_grow_pop.png",f)


# Compute multiple simulations of the evolution of generations
# Plot with error bands

f = GeneString.plot_multiple_epochs(pconst,sims_const)
save("monte_carlo_simulation_verify_hardy_weinberg_multi_simm_const_pop.png",f)

f = GeneString.plot_multiple_epochs(pgrow,sims_grow)
save("monte_carlo_simulation_verify_hardy_weinberg_multi_sim_grow_pop.png",f)





#= Fitness allele prediction 
# not in use 
# Compute the allele frequency
# Params is an array such that 
# p,rf_AA,rf_Aa,rf_aa = params
function compute_theoretical_allele_fitness_frequency(params)
    q = 1 - params[:p]
    # Calculate the mean fitness of the population
    w_bar = params[:p]^2 * params[:rf_AA] + 2 * params[:p] * q * params[:rf_Aa] + q^2 * params[:rf_aa]
    # Calculate the new allele frequency
    (params[:p]^2 * params[:rf_AA] + params[:p] * q * params[:rf_Aa]) / w_bar
end
=#





params = (
    genotypes = ["AA", "Aa", "aa"],
    num_sims = 10,
    num_generations = 100, # Number of generations
    population_size = 100000, # Size of population
    p = π/8, # Frequency of allele of interest "A"
    rf_AA = 1.0, # Relative fitness of AA
    rf_Aa = 1.0, # Relative fitness of Aa or aA
    rf_aa = 1.0) # Relative fitness of aa 

evolved_gen = evolve_generations(params)
f = plot_single_evolution_sample_path(params,evolved_gen)
save("monte_carlo_simulation_verify_hardy_weinberg_single_sim_irrational_frequ.png",f)


sims = simulate_multiple_evolutions(params)
f = plot_multiple_epochs(params,sims)
save("monte_carlo_simulation_verify_hardy_weinberg_multi_sim_irrational_frequ.png",f)




params = (
    genotypes = ["AA", "Aa", "aa"],
    num_sims = 10,
    num_generations = 10, # Number of generations
    population_size = 2, # Size of population
    p = π/8, # Frequency of allele of interest "A"
    rf_AA = 1.0, # Relative fitness of AA
    rf_Aa = 1.0, # Relative fitness of Aa or aA
    rf_aa = 1.0) # Relative fitness of aa 


    # Simulation initial population
    generations = Vector{String}[]
    @chain params begin
        _[:genotypes],compute_genotype_HW_frequency(_[:p]),_[:population_size]
        sample_initial_population(_...)
        push!(generations,_)
    end
    # Get state from population
    stats_per_generation =  Dict{String, Dict{String}}[]
    push!(stats_per_generation,compute_genotype_counts_frequency(generations[1]))
    

    # Evolve generations
    for gen in Base.OneTo(params[:num_generations])
        @chain gen begin
            next_generation(generations[_])
            @aside @chain _ begin
                push!(generations,_)
            end
            push!(stats_per_generation,compute_genotype_counts_frequency(_))
        end
    end

    Dict(
        "generations" => generations, 
        "stats" => stats_per_generation)





        fields_ha =     

        eyeball_est_mn_count = [250,600,120,250,1000]

        sum(eyeball_est_mn_count ./ fields_ha)*(1/5)


        sum(eyeball_est_mn_count) / sum(fields_ha)