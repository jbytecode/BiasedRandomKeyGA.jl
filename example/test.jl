# Traveling Salesman Example 


using BiasedRandomKeyGA


distances = [
    0.0 10.0 15.0 20.0 17.0 40;
    10.0 0.0 35.0 25.0 28.0 30;
    15.0 35.0 0.0 30.0 20.0 25;
    20.0 25.0 30.0 0.0 22.0 28;
    17.0 28.0 20.0 22.0 0.0 26;
    40.0 30.0 25.0 28.0 26.0 0.0
]

function costfunction(genes)
    # Encoder is sortperm.
    decodedval = sortperm(genes)
    totalcost = 0.0
    p = length(decodedval)
    for i in 1:(p-1)
        totalcost += distances[decodedval[i], decodedval[i+1]]
    end
    totalcost += distances[decodedval[p], decodedval[1]]
    return totalcost
end


ga = BRKGA(
    100,   # population size
    6,    # chromosome size
    [
        make_uniform_crossover(0.7),
        make_pathrelinking_crossover(0.1)
    ], # Multiple crossover functions
    [0.7, 0.3], # Crossover function probabilities
    10,    # number of elites
    10,    # number of mutants
    costfunction # cost function
)

function test()
    pop = create_population(ga)

    pop = generations(ga, pop, 1000) # Run for 1000 generations

    best = argmin(c -> c.cost, pop)

    println("Best cost found: ", best.cost)
    println("Best solution (decoded): ", sortperm(best.genes))
end 
