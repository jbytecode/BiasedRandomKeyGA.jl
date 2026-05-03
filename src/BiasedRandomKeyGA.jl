module BiasedRandomKeyGA

export BRKGA, Chromosome, Population
export create_population, generate_mutants
export generation, generations
export evaluate!
export selectacrossoverfunction
export make_uniform_crossover
export make_pathrelinking_crossover
export make_blx_crossover
export make_sbx_crossover
export createsimplega


mutable struct Chromosome
    genes::Vector{Float64}
    cost::Float64
end




"""
    struct BRKGA

    The `BRKGA` struct represents a Biased Random-Key Genetic Algorithm (BRKGA) instance. 
    It encapsulates all the parameters and configurations necessary to run the algorithm, 
    including population size, chromosome size, crossover functions, probabilities for selecting 
    crossover functions, number of elite chromosomes, number of mutants, and the cost function used
    to evaluate the fitness of chromosomes.    

# Fields

- `population_size::Int`: The number of chromosomes in the population.
- `chromosome_size::Int`: The number of genes in each chromosome.
- `crossoverfunctions::Vector{Function}`: A vector of crossover functions that can be used to 
   generate offspring from parent chromosomes.
- `crossoveropprobs::Vector{Float64}`: A vector of probabilities corresponding to each crossover 
   function, indicating the likelihood of selecting each function during the crossover process.
- `numelites::Int`: The number of elite chromosomes that are preserved in each generation without modification.
- `nummutants::Int`: The number of mutant chromosomes that are generated randomly in each generation to 
   maintain diversity in the population.
- `costfn::Function`: A function that takes a chromosome's genes as input and returns a cost value, 
   which is used to evaluate the quality of the chromosome in the context of the optimization problem being solved.
   In a BRKGA search, genes are passed directly to the cost function, which is responsible for interpreting 
   the genes and calculating the corresponding cost based on the specific problem domain. For example, 
   in a permutation-based problem, the cost function would need to decode the random keys into a permutation, possibly
   using the sortperm function. 
"""
struct BRKGA
    population_size::Int
    chromosome_size::Int
    crossoverfunctions::Vector{Function}
    crossoveropprobs::Vector{Float64}
    numelites::Int
    nummutants::Int
    costfn::Function
end


# Type alias for a population of chromosomes
const Population = Vector{Chromosome}


# Constructor for Chromosome that initializes genes with random values and cost with Inf
Chromosome(n::Int) = Chromosome(rand(n), Inf)



"""
    createsimplega(costfn::Function, chromosomesize::Int)::BRKGA

    Creates a simple BRKGA instance with default parameters for population size, crossover functions,
    crossover probabilities, number of elites, and number of mutants. The user can specify the cost 
    function and chromosome size, while the other parameters are set to commonly used defaults for a 
    basic BRKGA implementation.

# Arguments

- `costfn::Function`: A function that takes a chromosome's genes as input and returns a cost value, 
   which is used to evaluate the quality of the chromosome in the context of the optimization problem being solved.
- `chromosomesize::Int`: The number of genes in each chromosome, which determines the dimensionality of the 
   olution space being explored by the BRKGA.

# Notes

The default parameters used in this function are:
- `population_size`: 100
- `crossoverfunctions`: A vector containing a uniform crossover function with alpha = 0.7 
   and a path-relinking crossover function with sigma = 0.1.
- `crossoveropprobs`: A vector of probabilities corresponding to each crossover function, 
   indicating the likelihood of selecting each function during the crossover process. 
   In this case, the uniform crossover function has an 80% chance of being selected, while the
   path-relinking crossover function has a 20% chance of being selected.
- `numelites`: 10
- `nummutants`: 10
"""
function createsimplega(costfn::Function, chromosomesize::Int)::BRKGA
    crossoverfunctions = [make_uniform_crossover(0.7), make_pathrelinking_crossover(0.1)]
    crossoveropprobs = [0.8, 0.2]
    return BRKGA(100, chromosomesize, crossoverfunctions, crossoveropprobs, 10, 10, costfn)
end


"""

    make_blx_crossover(alpha::Float64)::Function

    Creates a BLX-alpha crossover function with the given alpha parameter. 
    The returned function takes two parent chromosomes and generates a child chromosome 
    by creating genes that are randomly selected from a range defined by the corresponding 
    genes of the two parents, expanded by a factor of alpha. This allows for exploration 
    of the solution space around the parents while maintaining some bias towards their gene values.

# Arguments

- `alpha::Float64`: A value between 0 and 1 that determines the extent of the range from which 
   the child genes are selected. A higher alpha allows for a wider range of possible gene values, 
   while a lower alpha restricts the child genes to be closer to the parent genes.

# Returns

- `Function`: A crossover function that can be used in the BRKGA framework to generate offspring 
   based on the BLX-alpha crossover method between two parent chromosomes.
"""
function make_blx_crossover(alpha=0.5)
    return (parent1, parent2) -> begin
        n = length(parent1)
        child = similar(parent1)
        for i in 1:n
            x1, x2 = parent1[i], parent2[i]
            d = abs(x1 - x2)
            lower = min(x1, x2) - alpha * d
            upper = max(x1, x2) + alpha * d
            # Clamp the child gene to be within [0, 1]
            child[i] = clamp(lower + rand() * (upper - lower), 0.0, 1.0)
        end
        return child
    end
end



"""
    make_sbx_crossover(nc::Float64)::Function

    Creates a Simulated Binary Crossover (SBX) function with the given distribution index (nc). 
    The returned function takes two parent chromosomes and generates a child chromosome by simulating 
    the binary crossover process in a continuous search space. The distribution index controls the 
    spread of the offspring around the parents, with a higher nc resulting in offspring that are 
    closer to the parents, and a lower nc allowing for more diverse offspring.

# Arguments

- `nc::Float64`: The distribution index for the SBX crossover, which controls the spread of the 
   offspring around the parents. A higher nc results in offspring that are closer to the parents, 
   while a lower nc allows for more diverse offspring.

# Returns

- `Function`: A crossover function that can be used in the BRKGA framework to generate offspring 
   based on the Simulated Binary Crossover method between two parent chromosomes.
"""
function make_sbx_crossover(nc = 20.0)
    return (p1, p2) -> begin
        n = length(p1)
        c1 = similar(p1)
        c2 = similar(p2)
        
        for i in 1:n
            if rand() <= 0.5 && abs(p1[i] - p2[i]) > 1e-9
                u = rand()
                if u <= 0.5
                    betaq = (2u)^(1/(nc + 1))
                else
                    betaq = (1/(2*(1-u)))^(1/(nc + 1))
                end
                
                c1[i] = clamp(0.5 * ((1 + betaq)*p1[i] + (1 - betaq)*p2[i]), 0.0, 1.0)
                c2[i] = clamp(0.5 * ((1 - betaq)*p1[i] + (1 + betaq)*p2[i]), 0.0, 1.0)
            else
                # If the parents are very close or by random chance, just copy the genes
                c1[i], c2[i] = p1[i], p2[i]
            end
        end
        return rand() < 0.5 ? c1 : c2
    end
end

"""
    make_uniform_crossover(alpha::Float64)::Function

Creates a uniform crossover function with the given alpha parameter.
This alpha parameter controls the bias towards the elite chromosome during crossover. 

# Arguments
- `alpha::Float64`: A value between 0 and 1 that determines the probability of inheriting a 
   gene from the elite chromosome. A higher alpha means a stronger bias towards the elite 
   chromosome, while a lower alpha allows for more diversity in the offspring.

# Returns
- `Function`: A crossover function that can be used in the BRKGA framework to generate offspring 
   based on a uniform crossover between an elite chromosome and another chromosome.
"""
function make_uniform_crossover(alpha::Float64)::Function
    return (ga::BRKGA, elitistc::Chromosome, c::Chromosome) -> begin
        n = length(elitistc.genes)
        child_genes = [rand() < alpha ? elitistc.genes[i] : c.genes[i] for i in 1:n]
        return Chromosome(child_genes, Inf)
    end
end





"""
    make_pathrelinking_crossover(sigma::Float64)::Function

Creates a path-relinking crossover function with the given sigma parameter.
The returned function takes a BRKGA instance, an elite chromosome, and another chromosome, 
and generates a new chromosome by exploring the path between the two parent chromosomes in the 
solution space. The offspring is generated by moving from the other chromosome towards the 
elite chromosome in steps defined by sigma, evaluating the cost at each step, and selecting 
the best solution found along the path.

# Arguments

- `sigma::Float64`: A parameter that controls the step size along the path between the two parent 
  chromosomes. A smaller sigma results in more fine-grained exploration, while a larger sigma allows 
  for more aggressive moves towards the elite chromosome.

# Returns

- `Function`: A crossover function that can be used in the BRKGA framework to generate offspring 
   based on path-relinking between an elite chromosome and another chromosome.

# References

- Noronha, Thiago F., and Celso C. Ribeiro. "Biased random-key genetic algorithms: 
  A tutorial with applications." Proceedings of the 2024 8th International Conference 
  on Intelligent Systems, Metaheuristics & Swarm Intelligence. 2024.
"""
function make_pathrelinking_crossover(sigma::Float64)::Function
    return (ga::BRKGA, elite::Chromosome, other::Chromosome) -> begin
        kmax = Int(floor(1 / sigma))
        bestgenes = Array{Float64}(undef, ga.chromosome_size)
        offgenes = Array{Float64}(undef, ga.chromosome_size)
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

# Function to generate mutant chromosomes
function generate_mutants(ga::BRKGA)::Vector{Chromosome}
    mutants = Vector{Chromosome}(undef, ga.nummutants)
    for i in 1:ga.nummutants
        mutants[i] = Chromosome(ga.chromosome_size)
    end
    return mutants
end

# Function to create an initial population of chromosomes
function create_population(ga::BRKGA)::Population
    population = Vector{Chromosome}(undef, ga.population_size)
    for i in 1:ga.population_size
        population[i] = Chromosome(ga.chromosome_size)
    end
    return population
end

# Evaluate the cost of a single chromosome
function evaluate!(ga::BRKGA, chrom::Chromosome)::Nothing
    chrom.cost = ga.costfn(chrom.genes)
    return nothing
end

# Evaluate the cost of each chromosome in the population
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
    sort!(population, by=c -> c.cost)

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
