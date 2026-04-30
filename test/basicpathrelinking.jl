@testset "Basic Path Relinking Tests" verbose = true begin

    @testset "Basic Permutation with Path Relinking" verbose = true begin

        function costfn(genes)
            decodedval = sortperm(genes)
            expected = collect(1:length(genes))     # 1, 2, ..., n
            return sum(abs.(decodedval .- expected))
        end


        # Configuration.
        ga = BRKGA(
            100,   # population size
            10,    # chromosome size
            50,    # generations
            [
                make_uniform_crossover(0.7),
                make_pathrelinking_crossover(0.1)
            ], # Multiple crossover functions
            10,    # number of elites
            20,    # number of mutants
            costfn # cost function
        )

        # Run the genetic algorithm.
        # Result is the final population with random keys.
        result = solve(ga)

        # Test the best solution 
        best = result[1]

        # Test if the best solution has zero cost 
        @test iszero(best.cost)

        # Test if the decoded solution is the expected permutation
        decoded_solution = sortperm(best.genes)
        expected_solution = collect(1:length(best.genes))
        @test decoded_solution == expected_solution # 1, 2, ..., n

    end


    @testset "1:50 Permutation" verbose = true begin

        function costfn(genes)
            decodedval = sortperm(genes)
            expected = collect(1:length(genes))     # 1, 2, ..., n
            return sum(abs.(decodedval .- expected))
        end


        # Configuration.
        ga = BRKGA(
            100,   # population size
            50,    # chromosome size
            1,    # generations
            [
                make_uniform_crossover(0.7),
                make_pathrelinking_crossover(0.1)
            ], # Multiple crossover functions
            10,    # number of elites
            10,    # number of mutants
            costfn # cost function
        )

        population = create_population(ga)
        iter = 0
        best = population[1]
        while true
            population = generation(ga, population)
            best = argmin(c -> c.cost, population)
            if iszero(best.cost)
                @info "Optimal solution found in iteration $iter"
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
            50,    # generations
            [
                make_uniform_crossover(0.7),
                make_pathrelinking_crossover(0.1)
            ], # Multiple crossover functions
            10,    # number of elites
            20,    # number of mutants
            costfn # cost function
        )

        result = solve(ga)
        best = result[1]
        bestperm = sortperm(best.genes)

        @test best.cost == 85.0

        # (1, 2, 3, 4), (4, 1, 2, 3), (3, 4, 1, 2) ... are all okay 
        # So I don't want to test for a specific permutation. 

    end
end