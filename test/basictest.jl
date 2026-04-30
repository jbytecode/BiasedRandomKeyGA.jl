@testset "Basic Tests" verbose = true begin

    @testset "Basic Permutation" verbose = true begin

        function costfn(genes)
            decodedval = sortperm(genes)
            expected = collect(1:length(genes))     # 1, 2, ..., n
            return sum(abs.(decodedval .- expected))
        end


        # Configuration.
        ga = BRKGA(
            100,   # population size
            10,    # chromosome size
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            10,    # number of elites
            20,    # number of mutants
            costfn # cost function
        )

       
        # Run the genetic algorithm.
        # Result is the final population with random keys.
        pop = create_population(ga)
        iter = 0
        best = pop[1]
        while true 
            pop = generation(ga, pop)
            evaluate!(ga, pop) # Ensure costs are updated after generation
            best = argmin(c -> c.cost, pop)
            if iszero(best.cost)
                break   
            end 
            iter += 1
            if iter > 10000
                @warn "Failed to find the optimal solution within 10000 iterations."
                break
            end 
        end 
        @test iszero(best.cost)
    end


    @testset "Big (50) Permutation" verbose = true begin

        function costfn(genes)
            decodedval = sortperm(genes)
            expected = collect(1:length(genes))     # 1, 2, ..., n
            return sum(abs.(decodedval .- expected))
        end


        # Configuration.
        ga = BRKGA(
            100,   # population size
            50,    # chromosome size
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            20,    # number of elites
            20,    # number of mutants
            costfn # cost function
        )

        population = create_population(ga)
        iter = 0
        best = population[1]
        while true 
            population = generation(ga, population)
            evaluate!(ga, population) # Ensure costs are updated after generation
            best = argmin(c -> c.cost, population)
            if iszero(best.cost)
                break
            end
            iter += 1
            if iter > 10000
                @warn "Failed to find the optimal solution within 10000 iterations."
                break
            end
        end 
        @test iszero(best.cost)

    end




    @testset "Basic Traveling Salesman" verbose = true begin

        # Distance matrix for 4 cities 
        distances = [0.0 10.0 15.0 20.0;
            10.0 0.0 25.0 45.0;
            15.0 25.0 0.0 30.0;
            20.0 45.0 30.0 0.0]

        # Cost function for TSP
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

        # GA Configuration for TSP
        ga = BRKGA(
            100,   # population size
            4,     # chromosome size (4 cities)
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            10,    # number of elites
            20,    # number of mutants
            costfn # cost function
        )
        
        pop = create_population(ga)
        iter = 0
        best = pop[1]
        while true
            pop = generation(ga, pop)
            evaluate!(ga, pop) # Ensure costs are updated after generation
            best = argmin(c -> c.cost, pop)
            if best.cost == 85.0
                break
            end
            iter += 1
            if iter > 10000
                @warn "Failed to find the optimal solution within 10000 iterations."
                @warn "The best solution found has cost $(best.cost)"
                @warn "Decoded solution: $(sortperm(best.genes))"
                break
            end
        end
        @test best.cost == 85.0
        

        # (1, 2, 3, 4), (4, 1, 2, 3), (3, 4, 1, 2) ... are all okay 
        # So I don't want to test for a specific permutation. 

    end
end