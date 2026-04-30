module BiasedRandomKeyGA

export BRKGA, Chromosome, Population
export create_population, generate_mutants
export generation, generations
export evaluate!
export selectacrossoverfunction
export make_uniform_crossover
export make_pathrelinking_crossover


mutable struct Chromosome
    genes::Vector{Float64}
    cost::Float64
end 


struct BRKGA
    population_size::Int
    chromosome_size::Int
    crossoverfunctions::Vector{Function}
    crossoveropprobs::Vector{Float64}
    numelites::Int
    nummutants::Int
    costfn::Function
end

const Population = Vector{Chromosome}

Chromosome(n::Int) = Chromosome(rand(n), Inf)

function make_uniform_crossover(alpha::Float64)::Function
    return (ga::BRKGA, elitistc::Chromosome, c::Chromosome) -> begin
        n = length(elitistc.genes)
        child_genes = [rand() < alpha ? elitistc.genes[i] : c.genes[i] for i in 1:n]
        return Chromosome(child_genes, Inf)
    end
end

function make_pathrelinking_crossover(sigma::Float64)::Function
    return (ga::BRKGA, elite::Chromosome, other::Chromosome) -> begin
        kmax = Int(floor(1 / sigma))
        bestgenes = Array{Float64}(undef, ga.chromosome_size)
        offgenes  = Array{Float64}(undef, ga.chromosome_size)
        bestcost = Inf

        for k in 1:kmax
            offgenes .= other.genes .+ (k * sigma) * (elite.genes .- other.genes)
            cost = ga.costfn(offgenes)
            if cost < bestcost
                bestcost = cost
                bestgenes .= offgenes
            end
        end
        return Chromosome(bestgenes, bestcost)
    end
end

function generate_mutants(ga::BRKGA)::Vector{Chromosome}
    mutants = Vector{Chromosome}(undef, ga.nummutants)
    for i in 1:ga.nummutants
        mutants[i] = Chromosome(ga.chromosome_size)
    end
    return mutants
end 

function create_population(ga::BRKGA)::Population
    population = Vector{Chromosome}(undef, ga.population_size)
    for i in 1:ga.population_size
        population[i] = Chromosome(ga.chromosome_size)
    end
    return population
end

function evaluate!(ga::BRKGA, chrom::Chromosome)::Nothing
    chrom.cost = ga.costfn(chrom.genes)
    return nothing
end

function evaluate!(ga::BRKGA, population::Population)::Nothing
    for chrom in population
        evaluate!(ga, chrom)
    end
    return nothing
end

function selectacrossoverfunction(ga::BRKGA)::Function
    weights = ga.crossoveropprobs
    cum_weights = cumsum(weights)
    r = rand()
    for (i, cw) in enumerate(cum_weights)
        if r < cw
            return ga.crossoverfunctions[i]
        end
    end
    return ga.crossoverfunctions[end] # Fallback
end

function generation(ga::BRKGA, population::Population)::Population
    # Crossover function
    crossfn = selectacrossoverfunction(ga)

    # Evaluate the cost of each chromosome
    evaluate!(ga, population)

    # Sort the population by cost
    sort!(population, by = c -> c.cost)

    # Select elites
    elites = population[1:ga.numelites]

    # Generate mutants
    mutants = generate_mutants(ga)

    new_population = Vector{Chromosome}(undef, ga.population_size)

    copyto!(new_population, 1, elites, 1, ga.numelites)
    copyto!(new_population, ga.numelites + 1, mutants, 1, ga.nummutants)

    nextindex = ga.numelites + ga.nummutants + 1
    for i in nextindex:ga.population_size
        elitistc = rand(elites)
        c = rand(population[ga.numelites+1:end]) # Select from non-elites
        new_population[i] = crossfn(ga, elitistc, c)
    end

    return new_population
end 

function generations(ga::BRKGA, population::Population, numgens::Int)::Population
    for i in 1:numgens
        population = generation(ga, population)
    end
    return population
end 



end # module BiasedRandomKeyGA
