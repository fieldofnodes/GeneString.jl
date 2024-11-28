module GeneString
    
    using Base.Iterators
    using Chain
    using StatsBase
    using CairoMakie
    using Random
    using Statistics
    
export
    AbstractOffspring,
    RandomNumberOffspring,
    FixedNumberOffspring,
    total_genome_single_gene,
    compute_genotype_HW_frequency,
    compute_genotype_counts_frequency,
    produce_punnet_square,
    convert_casing_upper_then_lower,
    produce_1_offspring,
    isHomozygous,
    isHeterozygous,
    isHomozygousDominant,
    isHomozygousRecessive,
    detectZygosity,
    split_genotype_into_pairs_iterator,
    sample_initial_population,
    number_offspring,
    next_generation,
    evolve_generations,
    evolve_stats,
    simulate_multiple_epochs,
    plot_single_epoch_sample_path,
    plot_multiple_epochs,
    std_error,
    AbstractInsect,
    GenericInsec

    std_error(x) = std(x)/length(x)


    abstract type AbstractInsect end
    struct GenericInsect <: AbstractInsect  
        id::Int
        sex::String
        gene::String
    end



    abstract type AbstractOffspring end

    """
            struct RandomNumberOffspring

        A structure to represent a range of random numbers for offspring generation.

        # Fields
        - `min::Number`: The minimum value of the range.
        - `max::Number`: The maximum value of the range.

        # Example
        ```julia
        expected_offspring_number = RandomNumberOffspring(1, 100)
        # Output: RandomNumberOffspring(1, 100)
        ```
    """
    struct RandomNumberOffspring <: AbstractOffspring
        min::Int
        max::Int    
    end



    """
            struct FixedNumberOffspring <: AbstractOffspring

        A structure to represent a fixed number for offspring generation.

        # Fields
        - `num::Number`: The number of offspring per reproduction

        # Example
        ```julia
        expected_offspring_number = FixedNumberOffspring(100)
        # Output: RandomNumberOffspring(100)
        ```
        
    """
    struct FixedNumberOffspring <: AbstractOffspring
        num::Int
    end


    """
            total_genome_single_gene(gene::AbstractString) -> Vector{String}

        Given a gene string `gene` with exactly two alleles, this function generates all 
        possible valid genome combinations (upper- and lowercase combinations of the alleles). 
        The resulting genome consists of the dominant-recessive combinations: `["AA", "Aa", "aa"]`, 
        with invalid combinations like `["aA"]` excluded.

        # Arguments:
        - `gene`: A string of length 2, where each character represents an allele (e.g., "Aa").

        # Returns:
        - A vector of strings containing all valid combinations of the genome.

        # Example:
        ```julia
        total_genome_single_gene("Aa")
        # Output: ["AA", "Aa", "aa"]
        ```
    """
    function total_genome_single_gene(gene::String)::Vector{String}
        @assert length(gene) == 2 "Function: total_genome_single_gene only suited for single gene"
       
        genome = @chain gene begin
            collect
            [lowercase(_[1]),uppercase(_[1])]
            Iterators.product(_,_)
            collect 
            join.(_) 
            reduce(vcat,_)
            sort
        end

        @chain genome begin
            collect.(_)  
            [!(islowercase(i[1]) && isuppercase(i[2])) for i in _] 
            genome[_]
        end
    end


    """
            compute_genotype_HW_frequency(dominant_allele_frequency::Float64) -> Vector{Float64}

        This function computes the theoretical genotype frequencies under Hardy-Weinberg 
        equilibrium for a population. Given the dominant allele frequency, it returns the 
        frequencies of homozygous dominant (AA), heterozygous (Aa), and homozygous recessive (aa) 
        genotypes.

        # Arguments:
        - `dominant_allele_frequency`: The frequency of the dominant allele in the population (a float between 0 and 1).

        # Returns:
        - A vector of floats representing the genotype frequencies in the following order: 
        - Homozygous dominant (AA)
        - Heterozygous (Aa)
        - Homozygous recessive (aa)

        # Example:
        ```julia
        compute_genotype_HW_frequency(0.7)
        # Output: [0.49, 0.42, 0.09]
        ```
    """ 
    function compute_genotype_HW_frequency(dominant_allele_frequency)
        recessive_allele_frequency = 1 - dominant_allele_frequency
        weight_homozygous_dominant = dominant_allele_frequency^2
        weight_heterozygous = 2*dominant_allele_frequency*recessive_allele_frequency
        weight_homozygous_recessive = recessive_allele_frequency^2
        [
            weight_homozygous_dominant, 
            weight_heterozygous, 
            weight_homozygous_recessive
        ]
    end


    """
            compute_genotype_counts_frequency(generation::Vector{String}) -> Dict{String, Dict{String, Float64}}

        This function computes both the genotype counts and the corresponding frequencies for a population 
        across a generation of genotypes. It compares the observed genotypes with the expected set based on a single gene 
        and returns both the counts and frequencies.

        # Arguments:
        - `generation`: A vector of strings where each string represents a genotype from a generation.

        # Returns:
        - A dictionary containing two sub-dictionaries:
        - `"counts"`: The counts of each observed genotype in the population.
        - `"frequencies"`: The relative frequencies of each genotype as a proportion of the total population.

        # Example:
        ```julia
        compute_genotype_counts_frequency(["AA", "Aa", "aa", "Aa", "AA"])
        # Output: Dict("counts" => Dict("AA" => 2, "Aa" => 2, "aa" => 1), "frequencies" => Dict("AA" => 0.4, "Aa" => 0.4, "aa" => 0.2))
        ```
    """
    function compute_genotype_counts_frequency(generation::Vector{String})::Union{Dict{String, Dict{String}},Dict{String, Int64}}
        if length(generation) == 0 
            @warn "No generation. Returning 0 for count and frequency"
            return Dict("counts" => 0, "frequencies" => 0)
        end
        expected_genome_set = total_genome_single_gene(generation[1])
        observed_genome_set = unique(generation)
        missing_genes = setdiff(expected_genome_set,observed_genome_set)
        counts = [count(ind -> ind == i, generation) for i in observed_genome_set]
        allele_counts = Dict(observed_genome_set[i] => counts[i] for i in eachindex(observed_genome_set))
        if missing_genes |> length != 0 # Some genes are missing
            [allele_counts[i] = 0 for i in missing_genes]
        end

        freqs = Base.values(allele_counts) ./ length(generation)
        allele_freq = Dict(collect(Base.keys(allele_counts))[i] => freqs[i] for i in eachindex(expected_genome_set))
        Dict("counts" => allele_counts, "frequencies" => allele_freq)
    end



    """
            produce_punnet_square(parent1::String, parent2::String)

        Generate a Punnett square from two parent genotypes.

        # Arguments
        - `parent1::String`: The genotype of the first parent.
        - `parent2::String`: The genotype of the second parent.

        # Returns
        - `Array{Tuple{Char,Char},1}`: An array of tuples representing the possible genotypes of the offspring.

        # Example
        ```julia
        produce_punnet_square("Aa", "Bb")
        # Output: [('A', 'B'), ('A', 'b'), ('a', 'B'), ('a', 'b')]
        ```
    """    
    function produce_punnet_square(parent1::String,parent2::String)
        @chain parent1,parent2 begin
            collect.(_)
            Iterators.product(_...)
            collect
        end
    end



    """
            convert_casing_upper_then_lower(gene::Tuple{Char, Char})

        Convert the casing of a gene tuple where the first character is lowercase and the second is uppercase.

        # Arguments
        - `gene::Tuple{Char, Char}`: A tuple representing a gene with two characters.

        # Returns
        - `Tuple{Char, Char}`: A tuple with the first character converted to uppercase and the second to lowercase if the original tuple had the first character in lowercase and the second in uppercase. Otherwise, returns the original tuple.

        # Example
        ```julia
        convert_casing_upper_then_lower(('a', 'B'))
        # Output: ('A', 'b')

        convert_casing_upper_then_lower(('A', 'b'))
        # Output: ('A', 'b')
        ```
    """
    function convert_casing_upper_then_lower(gene::Tuple)
        @assert length(gene) == 2 "Only meant to be length 2 and not $length(gene)."
        islowercase(gene[1]) & isuppercase(gene[2]) ? 
                (uppercase(gene[1]),lowercase(gene[2])) : 
                gene
    end


    """
            produce_1_offspring(parent1::String, parent2::String)

        Generate one offspring genotype from two parent genotypes using a Punnett square.

        # Arguments
        - `parent1::String`: The genotype of the first parent.
        - `parent2::String`: The genotype of the second parent.

        # Returns
        - `String`: A string representing the genotype of one possible offspring.

        # Example
        ```julia
        produce_1_offspring("Aa", "Bb")
        # Possible output: "AB", "Ab", "aB", or "ab"
        ```
    """
    function produce_1_offspring(parent1::String,parent2::String)
        @chain parent1,parent2 begin
            produce_punnet_square(_...)  
            sample
            convert_casing_upper_then_lower
            join
        end
    end

    """
            isHomozygous(gene::String)

        Determine if a gene is homozygous, meaning both alleles are the same case (either both uppercase or both lowercase).

        # Arguments
        - `gene::String`: A string representing a gene with two characters.

        # Returns
        - `Bool`: `true` if the gene is homozygous, `false` otherwise.

        # Example
        ```julia
        isHomozygous("AA")
        # Output: true

        isHomozygous("Aa")
        # Output: false

        ```
    """
    function isHomozygous(gene) 
        @assert length(gene) == 2 "Homozygous/Heterozygous only considered for diploid (AA or BB or Cd etc), not for gene of length $(length(gene))"
        cgene = collect(gene)
        all(isuppercase.(cgene)) || all(islowercase.(cgene))
    end

    

    """
            isHeterozygous(gene::String)

        Determine if a gene is heterozygous, meaning it contains one uppercase and one lowercase allele.

        # Arguments
        - `gene::String`: A string representing a gene with two characters.

        # Returns
        - `Bool`: `true` if the gene is heterozygous, `false` otherwise.

        # Example
        ```julia
        isHeterozygous("Aa")
        # Output: true

        isHeterozygous("AA")
        # Output: false
        ```
    """
    function isHeterozygous(gene) 
        @assert length(gene) == 2 "Homozygous/Heterozygous only considered for diploid (AA or BB or Cd etc), not for gene of length $(length(gene))"
        cgene = collect(gene)
        
        hetero_bool = any(isuppercase.(cgene)) && any(islowercase.(cgene))
        @assert hetero_bool == !isHomozygous(gene) "Genes are meant to equal the negation of the other for Homozygous versus Heterozygous, they do not. Check fucntions: isHomozygous and isHeterozygous" 
        hetero_bool
    end




    """
            isHomozygousDominant(gene::String)

        Determine if a gene is homozygous dominant, meaning both alleles are uppercase.

        # Arguments
        - `gene::String`: A string representing a gene with two characters.

        # Returns
        - `Bool`: `true` if the gene is homozygous dominant, `false` otherwise.

        # Example
        ```julia
        isHomozygousDominant("AA")
        # Output: true

        isHomozygousDominant("Aa")
        # Output: false
        ```
    """
    function isHomozygousDominant(gene)
        @assert length(gene) == 2 "Homozygous/Heterozygous only considered for diploid (AA or BB or Cd etc), not for gene of length $(length(gene))"
        if !(isHomozygous(gene))
            return false
        end
        cgene = collect(gene)
        all(isuppercase.(cgene)) 
    end




    """
            isHomozygousRecessive(gene::String)

        Determine if a gene is homozygous recessive, meaning both alleles are lowercase.

        # Arguments
        - `gene::String`: A string representing a gene with two characters.

        # Returns
        - `Bool`: `true` if the gene is homozygous recessive, `false` otherwise.

        # Example
        ```julia
        isHomozygousRecessive("aa")
        # Output: true

        isHomozygousRecessive("Aa")
        # Output: false
        ```
    """
    function isHomozygousRecessive(gene)
        @assert length(gene) == 2 "Homozygous/Heterozygous only considered for diploid (AA or BB or Cd etc), not for gene of length $(length(gene))"
        if !(isHomozygous(gene))
            return false
        end
        cgene = collect(gene)
        hom_rec = all(islowercase.(cgene))
        @assert hom_rec == !isHomozygousDominant(gene) "Genes are meant to equal the negation of the other for Homozygous versus Heterozygous, they do not. Check fucntions: isHomozygous and isHeterozygous" 
        hom_rec
    end


    """
            detectZygosity(gene::String)

        Determine the zygosity of a gene, identifying it as homozygous dominant, homozygous recessive, or heterozygous.

        # Arguments
        - `gene::String`: A string representing a gene with two characters.

        # Returns
        - `String`: A string indicating the zygosity of the gene: "HomozygousDominant", "HeterozygousRecessive", or "Heterozygous".

        # Example
        ```julia
        detectZygosity("AA")
        # Output: "HomozygousDominant"

        detectZygosity("aa")
        # Output: "HeterozygousRecessive"

        detectZygosity("Aa")
        # Output: "Heterozygous"
        ```
    """
    function detectZygosity(gene)
        if isHomozygous(gene)
            if isHomozygousDominant(gene)
                return "HomozygousDominant"
            elseif isHomozygousRecessive(gene)
                return "HeterozygousRecessive"
            else 
                error("The case of being homozygous and not being dominant or recessive should not happen, check functions: isHomozygousRecessive, isHomozygousRecessive")
            end
        elseif isHeterozygous(gene) 
            return "Heterozygous"
        else
            error("The case of not being homozygous not heterozygous should not happen, check functions: isHomozygousRecessive, isHomozygousRecessive")
        end
    end

    """
        split_genotype_into_pairs_iterator(genotype::String)

    Split a genotype string into pairs of characters using an iterator.

    # Arguments
    - `genotype::String`: A string representing a genotype.

    # Returns
    - `Base.Iterators.Partition`: An iterator that yields pairs of characters from the genotype string.

    # Example
    ```julia
    split_genotype_into_pairs_iterator("AABB")
    # Output: Base.Iterators.Partition{String, Tuple{SubString{String}, SubString{String}}}("AABB", 2)
    ```
    """
    function split_genotype_into_pairs_iterator(genotype)
        Iterators.partition(genotype,2)
    end
    
    """
            sample_initial_population(genotypes::Vector{String}, weights::Vector{Float64}, population_size::Int)

        Generate an initial population sample from given genotypes and their corresponding weights.

        # Arguments
        - `genotypes::Vector{String}`: A vector of genotype strings.
        - `weights::Vector{Float64}`: A vector of weights corresponding to the genotypes.
        - `population_size::Int`: The desired size of the population sample.

        # Returns
        - `Vector{String}`: A vector of genotype strings representing the sampled population.

        # Example
        ```julia
        genotypes = ["AA", "Aa", "aa"]
        weights = [0.5, 0.3, 0.2]
        population_size = 10
        sample_initial_population(genotypes, weights, population_size)
        # Possible output: ["AA", "Aa", "AA", "aa", "AA", "Aa", "AA", "AA", "Aa", "aa"]
        ```
    """
    function sample_initial_population(genotypes,weights,population_size)
        sample(genotypes,Weights(weights), population_size)
    end

    





    """
            number_offspring(offspring_struct::AbstractOffspring)

        Abstract function to determine the number of offspring. This function should be implemented by specific offspring types.
    """
    function number_offspring(offspring_struct::AbstractOffspring) end


    """
            number_offspring(offspring_struct::RandomNumberOffspring) -> Int

        Returns a random number of offspring within the specified range of the `RandomNumberOffspring` struct.

        # Arguments
        - `offspring_struct::RandomNumberOffspring`: A struct containing the minimum and maximum number of offspring.

        # Returns
        - `Int`: A random number of offspring between `offspring_struct.min` and `offspring_struct.max`.
    """
    function number_offspring(offspring_struct::RandomNumberOffspring) 
        rand(offspring_struct.min:offspring_struct.max)
    end


    """
            number_offspring(offspring_struct::FixedNumberOffspring) -> Int

        Returns a fixed number of offspring as specified in the `FixedNumberOffspring` struct.

        # Arguments
        - `offspring_struct::FixedNumberOffspring`: A struct containing the fixed number of offspring.

        # Returns
        - `Int`: The fixed number of offspring specified by `offspring_struct.num`.
    """
    function number_offspring(offspring_struct::FixedNumberOffspring) 
        offspring_struct.num
    end
        


    """
    next_generation(ran_off::RandomNumberOffspring, population::Vector{String})

    Generate the next generation of a population by sampling pairs of parents and producing a random number of offspring within a specified range.

    # Arguments
    - `ran_off::RandomNumberOffspring`: An instance of `RandomNumberOffspring` specifying the range for the number of offspring.
    - `population::Vector{String}`: A vector of genotype strings representing the current population.

    # Returns
    - `Vector{String}`: A vector of genotype strings representing the next generation.

    # Example
    ```julia
    ran_off = RandomNumberOffspring(1, 3)
    current_population = ["AA", "Aa", "aa", "Aa"]
    next_generation(ran_off, current_population)
    # Possible output: ["AA", "Aa", "aa", "AA", "Aa", "aa"]
    ```
    """
    function next_generation(offspring_struct::AbstractOffspring ,population)
        update_generation = String[]
        for _ in Base.eachindex(population)
            parents = sample(population,2)
            num_offspring = number_offspring(offspring_struct)
            for offspring in Base.OneTo(num_offspring)
                offspring = produce_1_offspring(parents...)
                push!(update_generation,offspring)
            end
        end
        update_generation 
    end


    

    
    """
            evolve_generations(
                offspring_struct::AbstractOffspring,
                number_generations::Int,
                generations::Vector{Vector{String}}
            ) -> Vector{Vector{String}}

        Simulates the evolution of generations based on the given offspring structure and number of generations.

        # Arguments
        - `offspring_struct::AbstractOffspring`: The structure defining the offspring generation rules.
        - `number_generations::Int`: The number of generations to evolve.
        - `generations::Vector{Vector{String}}`: A vector containing the initial population as the first element.

        # Returns
        - `Vector{Vector{String}}`: A vector where each element is a generation of the population.

        # Example
        ```julia
        initial_population = [["AA", "Aa", "aa", "Aa"]]
        offspring_struct = RandomNumberOffspring(1, 3)
        evolved_generations = evolve_generations(offspring_struct, 5, initial_population)
        println(evolved_generations)

        ```
    """
    function evolve_generations(
        offspring_struct::AbstractOffspring,
        number_generations,
        generations::Vector{Vector{String}})

        @assert length(generations) == 1 "Variable: generations is the initial population only and should be length 1, not $(length(generations))"
        for gen in Base.OneTo(number_generations)
            @chain gen begin
                next_generation(offspring_struct,generations[_])
                push!(generations,_)
            end
        end
        generations
        
    end



    """
            evolve_stats(epoch::Vector{Dict{String, Any}}) -> Dict{String, Any}

        Computes and returns statistics for each generation in the given epoch.

        # Arguments
        - `epoch::Vector{Dict{String, Any}}`: A vector of dictionaries, where each dictionary represents the genotype frequencies for a generation.

        # Returns
        - `Dict{String, Any}`: A dictionary containing the following keys:
            - `"homozygous_dominant"`: A vector of frequencies for the homozygous dominant genotype across generations.
            - `"homozygous_recessive"`: A vector of frequencies for the homozygous recessive genotype across generations.
            - `"heterozygous"`: A vector of frequencies for the heterozygous genotype across generations.
            - `"generations_vec"`: A range object representing the generations.

        # Example
        ```julia
        epoch = [
            Dict("frequencies" => Dict("AA" => 0.5, "aa" => 0.2, "Aa" => 0.3)),
            Dict("frequencies" => Dict("AA" => 0.6, "aa" => 0.1, "Aa" => 0.3))
        ]
        stats = evolve_stats(epoch)
        println(stats)
        ```
    """
    function evolve_stats(epoch)
        num_gens = length(epoch)
        stats = compute_genotype_counts_frequency.(epoch)
        gene_letters = @chain stats begin
            _[1]["frequencies"]
            keys
            collect 
            collect.(_)
            reduce(vcat,_)
            unique
        end
        hd = uppercase.(gene_letters) |> join  
        hr = lowercase.(gene_letters) |> join  
        h = (uppercase(gene_letters[1]),lowercase(gene_letters[2])) |> join
             
        dict = Dict{String,Any}()
        dict["homozygous_dominant"] = [i["frequencies"][hd] for i in stats]
        dict["homozygous_recessive"] = [i["frequencies"][hr] for i in stats]
        dict["heterozygous"] = [i["frequencies"][h] for i in stats]
        dict["generations_vec"] = Base.OneTo(num_gens)
        dict
    end


    """
        simulate_multiple_epochs(
            offspring_struct::AbstractOffspring,
            number_generations::Int,
            number_sims::Int,
            genotypes::Vector{String},
            weights::Vector{Float64},
            population_size::Int
        ) -> Dict{Int, Dict{String, Any}}

        Simulates multiple epochs of population evolution based on the given parameters.

        # Arguments
        - `offspring_struct::AbstractOffspring`: The structure defining the offspring generation rules.
        - `number_generations::Int`: The number of generations to evolve in each simulation.
        - `number_sims::Int`: The number of simulations to run.
        - `genotypes::Vector{String}`: A vector of possible genotypes.
        - `weights::Vector{Float64}`: A vector of weights corresponding to the genotypes.
        - `population_size::Int`: The size of the initial population.

        # Returns
        - `Dict{Int, Dict{String, Any}}`: A dictionary where each key is a simulation index and each value is a dictionary containing:
            - `"generation"`: A vector of generations.
            - `"stats"`: A dictionary of statistics for each generation.

        # Example
        ```julia
        offspring_struct = RandomNumberOffspring(1, 3)
        number_generations = 5
        number_sims = 10
        genotypes = ["AA", "Aa", "aa"]
        weights = [0.5, 0.3, 0.2]
        population_size = 100

        simulation_history = simulate_multiple_epochs(
            offspring_struct, number_generations, number_sims, genotypes, weights, population_size
        )
        println(simulation_history)

        ```
    """
    function simulate_multiple_epochs(
    offspring_struct::AbstractOffspring,
    number_generations,
    number_sims,
    genotypes,
    weights,
    population_size)

    Simulation_history = Dict{Int, Any}()
    for sim in Base.OneTo(number_sims)
        generations = Vector{String}[]
        N₀ = sample_initial_population(genotypes,weights,population_size)
        push!(generations,N₀)
        epoch = evolve_generations(
            offspring_struct,number_generations,generations)
        stats = evolve_stats(epoch)
        Simulation_history[sim] = Dict("generation" => epoch,"stats" => stats)
    end
    @warn "2024-10-02: Outputting all generational data along with stats. Most likely remove specific generation results later and just leave stats"
    Simulation_history
    end





        

    """
        plot_single_epoch_sample_path(params::NamedTuple, evolved_gen::Vector{Vector{String}}) -> Figure

    Plots the evolution of genotype frequencies over generations for a single simulation sample path.

    # Arguments
    - `params::NamedTuple`: A named tuple containing the parameters for the simulation, including:
        - `:p`: The allele frequency.
        - `:population_size`: The size of the population.
        - `:num_generations`: The number of generations.
    - `evolved_gen::Vector{Vector{String}}`: A vector of vectors representing the evolved generations.

    # Returns
    - `Figure`: A figure object displaying the evolution of genotype frequencies.
        ```
    """
    function plot_single_epoch_sample_path(params::NamedTuple,evolved_gen::Vector{Vector{String}})
        off_type = params[:offspring] |> typeof |> string
        weights = compute_genotype_HW_frequency(params[:p])
        stats = evolve_stats(evolved_gen)
        f = Figure(size = (800,800),fontsize = 22)
        ax = Axis(f[1,1],aspect = 1, xlabel = "Generation", ylabel = "Frequency",title = "Monte Carlo simulation to verify Hardy-Weinberg principle\nN₀ = $(params[:population_size]), Sims = 1, Gens = $(params[:num_generations]), Offspring: $(off_type)")
        lines!(ax,stats["generations_vec"],stats["homozygous_dominant"], label = "Homozygous dominant",color = :blue)
        lines!(ax,stats["generations_vec"],stats["homozygous_recessive"], label = "Homozygous recessive", color = :green)
        lines!(ax,stats["generations_vec"],stats["heterozygous"], label = "Heterozygous", color = :purple)
        hlines!(ax,weights[1],label = "Theory",color = :red,linestyle = :dash)
        hlines!(ax,weights[2],color = :red,linestyle = :dash)
        hlines!(ax,weights[3],color = :red,linestyle = :dash)
        axislegend()
        ylims!(0,1)
        f
    end
   

    """
            plot_multiple_epochs(params::NamedTuple, sims::Dict{Int64, Any}) -> Figure

        Plots the evolution of genotype frequencies over multiple simulations.

        # Arguments
        - `params::NamedTuple`: A named tuple containing the parameters for the simulation, including:
            - `:p`: The allele frequency.
            - `:population_size`: The size of the population.
            - `:num_sims`: The number of simulations.
            - `:num_generations`: The number of generations.
            - `:offspring`: The offspring structure.
        - `sims::Vector{Dict{String, Any}}`: A vector of dictionaries, each representing the results of a single simulation.

        # Returns
        - `Figure`: A figure object displaying the evolution of genotype frequencies with mean and standard error bands.
        ```
    """
    function plot_multiple_epochs(params::NamedTuple,sims::Dict{Int64, Any})
    off_type = params[:offspring] |> typeof |> string
    weights = compute_genotype_HW_frequency(params[:p])

    num_gens = sims[1]["stats"]["generations_vec"]

    # Homozygous dominant
    homo_dom = [sims[i]["stats"]["homozygous_dominant"] for i in eachindex(sims) ] |> 
    x -> reduce(hcat,x) |>
    x -> x'
    mn_homozygous_dominant = mean(homo_dom,dims = 1) |> vec
    se_homozygous_dominant = [std_error(i) for i in eachcol(homo_dom)]
    lo_homozygous_dominant = mn_homozygous_dominant .- se_homozygous_dominant
    up_homozygous_dominant = mn_homozygous_dominant .+ se_homozygous_dominant

    # Homozygous recessive
    homo_rec = [sims[i]["stats"]["homozygous_recessive"] for i in eachindex(sims) ] |> 
    x -> reduce(hcat,x) |>
    x -> x'
    mn_homozygous_recessive = mean(homo_rec,dims = 1) |> vec
    se_homozygous_recessive = [std_error(i) for i in eachcol(homo_rec)]
    lo_homozygous_recessive = mn_homozygous_recessive .- se_homozygous_recessive
    up_homozygous_recessive = mn_homozygous_recessive .+ se_homozygous_recessive

    # Heterozygous
    hetero = [sims[i]["stats"]["heterozygous"] for i in eachindex(sims) ] |> 
    x -> reduce(hcat,x) |>
    x -> x'
    mn_heterozygous = mean(hetero,dims = 1) |> vec
    se_heterozygous = [std_error(i) for i in eachcol(hetero)]
    lo_heterozygous = mn_heterozygous .- se_heterozygous
    up_heterozygous = mn_heterozygous .+ se_heterozygous


    f = Figure(size = (800,800),fontsize = 22)
    ax = Axis(f[1,1],aspect = 1, xlabel = "Generation", 
        ylabel = "Frequency",
        title = "Monte Carlo simulation to verify Hardy-Weinberg principle\nN₀ = $(params[:population_size]), Sims = $(params[:num_sims]), Gens = $(params[:num_generations]), Offspring: $(off_type)")
    lines!(ax,num_gens,mn_homozygous_dominant,color = :blue,label = "Homozygous dominant") 
    lines!(ax,num_gens,mn_homozygous_recessive,color = :green,label = "Homozygous recessive") 
    lines!(ax,num_gens,mn_heterozygous,color = :purple,label = "Heterozygous") 
    band!(ax,num_gens,lo_homozygous_dominant,up_homozygous_dominant, color = (:blue,0.3), label = "Standard error")
    band!(ax,num_gens,lo_homozygous_recessive,up_homozygous_recessive, color = (:green,0.3))
    band!(ax,num_gens,lo_heterozygous,up_heterozygous, color = (:purple,0.3))
    hlines!(ax,weights[1],label = "Theory",color = :red,linestyle = :dash)
    hlines!(ax,weights[2],color = :red,linestyle = :dash)
    hlines!(ax,weights[3],color = :red,linestyle = :dash)
    axislegend()
    ylims!(0,1)
    f
end


function ascribe_sex_based_on_gene(population)
    gender_depend_gene = String[]
    for n in eachindex(population)
        if GeneString.isHomozygousDominant(population[n])
            push!(gender_depend_gene,"M")
        elseif GeneString.isHomozygousRecessive(population[n])
            push!(gender_depend_gene,rand(["M","F"]))
        elseif GeneString.isHeterozygous(population[n])
            push!(gender_depend_gene,rand(["M","F"]))
        else
            @error "Error no gene present matching your input"
        end
    end
    gender_depend_gene
end

function next_generation(offspring_struct::AbstractOffspring,population,sex)
    @assert length(population) == length(sex) "Population and sex must have same length"
    update_generation = String[]
    fem_inds = findall(x->x=="F",sex)
    mal_inds = findall(x->x=="M",sex)
    for mom_ind in fem_inds
        dad_ind = sample(mal_inds)
        num_offspring = number_offspring(offspring_struct)
        for offspring in Base.OneTo(num_offspring)
            offspring = produce_1_offspring(population[dad_ind],population[mom_ind])
            push!(update_generation,offspring)
        end
    end
    update_generation 
    
end

function indices_to_keep(population,sex)
    @assert length(population) == length(sex) "Population and sex must have same length"
    inds = eachindex(population)
    f_inds = findall(x->x=="F",sex)
    ox_gene_inds = findall(x->isHeterozygous(x)||isHomozygousDominant(x),population)
    removal = intersect(f_inds,ox_gene_inds)
    setdiff(inds,removal)
end



function count_unique(x)
    ux = unique(x)
    y = map(i -> count(==(i),x),ux)
    Dict(ux .=> y)
end





end






