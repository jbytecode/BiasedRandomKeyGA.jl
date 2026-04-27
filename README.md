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
            10,    # chromosome size
            50,    # generations
            0.7,   # alpha (for Uniform Crossover)
            10,    # number of elites
            20,    # number of mutants
            costfn # cost function
        )
```

The `solve` method calculates the steps:

```julia
solve(ga)
```

The `solve` method returns the sorted population. The best solution can be handled using 

```julia
best = result[1]
```

Since the chromosomes are encoded using real values, an `encoder` function can be 
used the transform original chromosome into a permutation vector:

```julia
decoded_solution = sortperm(best.genes)
```

A cost function would be in form of

```julia
function costfn(genes)
    decodedval = sortperm(genes)
    ...
    ...
    return totalcost
end
```

and the function always returns the cost. The first line transforms real numbers into a permutation
vector. A possible implementation of a cost function in a Traveling Salesman Problem seems like 

```julia
function costfn(genes)
    decodedval = sortperm(genes)
    totalcost = 0.0
    p = length(decodedval)
    for i in 1:(p-1)
        totalcost += distances[decodedval[i], decodedval[i+1]]
    end
    totalcost += distances[decodedval[p], decodedval[1]]
    return totalcost
end
```
