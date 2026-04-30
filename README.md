# BiasedRandomKeyGA.jl


Biased Random Key Genetic Algorithms with Julia


# Installation 

```julia
julia> ]
pkg> add https://github.com/jbytecode/BiasedRandomKeyGA.jl
```

# Example 

The implementation requires the parameters in a single struct `BRKGA`: 

```julia
        ga = BRKGA(
            100,   # population size
            50,    # chromosome size
            [
                make_uniform_crossover(0.7),
                make_pathrelinking_crossover(0.1)
            ], # Multiple crossover functions
            [0.7, 0.3], # Crossover function probabilities
            10,    # number of elites
            10,    # number of mutants
            costfn # cost function
        )
```

Creating an initial population:

```julia
        pop = create_population(ga)
```

Iterating the population: 

```julia
        pop = generation(ga, pop)
```

Iterating the population 100 times: 

```julia
        for i in 1:100
            pop = generation(ga, pop)
        end
```

Getting the final results:

```julia
        evaluate!(ga, pop) # Ensure costs are updated before finding the best solution
        best = argmin(c -> c.cost, pop)
```


A sample cost function:

```julia
function costfn(chromosome)
    # Decode first
    solution = sortperm(chromosome.genes)
    
    # Do calculations 
    # ....

    # Return the cost value. Objective is to minimize this value.
    return cost_value
end
```