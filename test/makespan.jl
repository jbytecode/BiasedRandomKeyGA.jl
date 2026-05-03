struct Process{T1 <: Real, T2 <: Real, T3 <: Real} 
    start::T1
    duration::T2
    finish::T3
end

function makespan(times::Matrix, permutation::Vector{Int})::Float64

    n, m = size(times)

    timetable = Matrix{Process}(undef, m, n)

    for machine_id = 1:m
        for task_id = 1:n
            current_task = permutation[task_id]
            if machine_id == 1
                if task_id == 1
                    start = 0
                else
                    start = timetable[machine_id, task_id-1].finish
                end
            else
                if task_id == 1
                    start = timetable[machine_id-1, task_id].finish
                else
                    start = max(
                        timetable[machine_id, task_id-1].finish,
                        timetable[machine_id-1, task_id].finish,
                    )
                end
            end
            duration = times[current_task, machine_id]
            finish = start + duration
            timetable[machine_id, task_id] = Process(start, duration, finish)
        end
    end

    return timetable[end, end].finish
end




@testset "Makespan test" verbose = true begin 

    @testset "Mini Makespan test - 1 (2-machines)" begin 
        
       times = Float64[
                3.2 4.2
                4.7 1.5
                2.2 5.0
                5.8 4.0
                3.1 2.8
            ]

        bestcost = 20.5

        function makespancostfn(genes::Vector{Float64})::Float64
            permutation = sortperm(genes)
            return makespan(times, permutation)
        end

        ga = BRKGA(
            100,   # population size
            5,    # chromosome size
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            20,    # number of elites
            20,    # number of mutants
            makespancostfn # cost function
        )

        pop = create_population(ga)
        maxiter = 10000
        iter = 0
        for i in 1:maxiter
            pop = generation(ga, pop)
            if iter % 100 == 0
                evaluate!(ga, pop) # Ensure costs are updated after generation
                best = argmin(c -> c.cost, pop)
                if best.cost == bestcost
                    break
                end
            end
            iter += 1
        end
        
        evaluate!(ga, pop) # Ensure costs are updated after generation
        
        best = argmin(c -> c.cost, pop)

        @test best.cost == bestcost
    end 




    @testset "Mini Makespan test - 2 (2-machines)" begin 
        
       times = Float64[
                4 7
                8 3
                5 8
                6 4
                8 5
                7 4
            ]

        bestcost = 41

        function makespancostfn(genes::Vector{Float64})::Float64
            permutation = sortperm(genes)
            return makespan(times, permutation)
        end

        ga = BRKGA(
            100,   # population size
            6,    # chromosome size
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            20,    # number of elites
            20,    # number of mutants
            makespancostfn # cost function
        )

        pop = create_population(ga)
        maxiter = 10000
        iter = 0
        for i in 1:maxiter
            pop = generation(ga, pop)
            if iter % 100 == 0
                evaluate!(ga, pop) # Ensure costs are updated after generation
                best = argmin(c -> c.cost, pop)
                if best.cost == bestcost
                    break
                end
            end
            iter += 1
        end
        
        evaluate!(ga, pop) # Ensure costs are updated after generation
        
        best = argmin(c -> c.cost, pop)

        @test best.cost == bestcost
    end




    @testset "Basic Makespan (3 machines)" begin 

        times = Float64[
                3 3 5
                8 4 8
                7 2 10
                5 1 7
                2 5 6
            ]

        bestcost = 42 

        function makespancostfn(genes::Vector{Float64})::Float64
            permutation = sortperm(genes)
            return makespan(times, permutation)
        end

        ga = BRKGA(
            100,   # population size
            5,    # chromosome size
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            20,    # number of elites
            20,    # number of mutants
            makespancostfn # cost function
        )

        pop = create_population(ga)
        maxiter = 10000
        iter = 0
        for i in 1:maxiter
            pop = generation(ga, pop)
            if iter % 100 == 0
                evaluate!(ga, pop) # Ensure costs are updated after generation
                best = argmin(c -> c.cost, pop)
                if best.cost == bestcost
                    break
                end
            end
            iter += 1
        end
        
        evaluate!(ga, pop) # Ensure costs are updated after generation
        
        best = argmin(c -> c.cost, pop)

        @test best.cost == bestcost
    end 




    @testset "Four machines" begin 

        times = Float64[
                3 3 5 2
                8 4 8 3
                7 2 10 4
                5 1 7 5
                2 5 6 6
            ]

        bestcost = 44

        function makespancostfn(genes::Vector{Float64})::Float64
            permutation = sortperm(genes)
            return makespan(times, permutation)
        end

        ga = BRKGA(
            100,   # population size
            5,    # chromosome size
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            20,    # number of elites
            20,    # number of mutants
            makespancostfn # cost function
        )

        pop = create_population(ga)
        maxiter = 10000
        iter = 0
        for i in 1:maxiter
            pop = generation(ga, pop)
            if iter % 100 == 0
                evaluate!(ga, pop) # Ensure costs are updated after generation
                best = argmin(c -> c.cost, pop)
                if best.cost == bestcost
                    break   
                end 
            end 
            iter += 1
        end 

        evaluate!(ga, pop) # Ensure costs are updated after generation

        best = argmin(c -> c.cost, pop)
    
        @test best.cost == bestcost

    end 



    @testset "Five machines" begin 

         times = Float64[
                7 5 2 3 9
                6 6 4 5 10
                5 4 5 6 8
                8 3 3 2 6
            ]

        bestcost = 51

        function makespancostfn(genes::Vector{Float64})::Float64
            permutation = sortperm(genes)
            return makespan(times, permutation)
        end

        ga = BRKGA(
            100,   # population size
            4,    # chromosome size
            [make_uniform_crossover(0.7)], # crossover function
            [1.0], # crossover function probabilities
            20,    # number of elites
            20,    # number of mutants
            makespancostfn # cost function
        )

        pop = create_population(ga)
        maxiter = 10000
        iter = 0
        for i in 1:maxiter
            pop = generation(ga, pop)
            if iter % 100 == 0
                evaluate!(ga, pop) # Ensure costs are updated after generation
                best = argmin(c -> c.cost, pop)
                if best.cost == bestcost
                    break   
                end 
            end 
            iter += 1
        end 

        evaluate!(ga, pop) # Ensure costs are updated after generation

        best = argmin(c -> c.cost, pop)
    
        @test best.cost == bestcost
    end 
end 