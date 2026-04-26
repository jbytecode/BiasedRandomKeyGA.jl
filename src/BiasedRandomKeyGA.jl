module BiasedRandomKeyGA

export BRKGA, Chromosome
export create_population, generation
export evaluate!, solve

mutable struct Chromosome
    genes::Vector{Float64}
    cost::Float64
end 


struct BRKGA
    population_size::Int
    chromosome_size::Int
    generations::Int
    alpha::Float64
    numelites::Int
    nummutants::Int
    costfn::Function
end

const Population = Vector{Chromosome}

Chromosome(n::Int) = Chromosome(rand(n), Inf)

function uniform_crossover(ga::BRKGA, elitistc::Chromosome, c::Chromosome)
    n = length(elitistc.genes)
    alpha = ga.alpha
    child_genes = [rand() < alpha ? elitistc.genes[i] : c.genes[i] for i in 1:n]
    return Chromosome(child_genes, Inf)
end 

function generate_mutants(ga::BRKGA)
    return [Chromosome(ga.chromosome_size) for _ in 1:ga.nummutants]
end 

function create_population(ga::BRKGA)
    return [Chromosome(ga.chromosome_size) for _ in 1:ga.population_size]
end 

function evaluate!(ga::BRKGA, population::Population)::Nothing
    for chrom in population
        chrom.cost = ga.costfn(chrom.genes)
    end
    return nothing
end

function generation(ga::BRKGA, population::Population)
    # Evaluate the cost of each chromosome
    evaluate!(ga, population)

    # Sort the population by cost
    sort!(population, by = c -> c.cost)

    # Select elites
    elites = population[1:ga.numelites]

    # Generate mutants
    mutants = generate_mutants(ga)

    # Create new population with elites and mutants
    new_population = vcat(elites, mutants)

    # Fill the rest of the population with offspring from crossover
    while length(new_population) < ga.population_size
        elitistc = rand(elites)
        c = rand(population[ga.numelites+1:end]) # Select from non-elites
        child = uniform_crossover(ga, elitistc, c)
        push!(new_population, child)
    end

    return new_population
end 



function solve(ga::BRKGA)
    population = create_population(ga)
    for _ in 1:ga.generations
        population = generation(ga, population)
    end
    evaluate!(ga, population) # Final evaluation to get the best solution
    population
end 







end # module BiasedRandomKeyGA
