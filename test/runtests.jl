using Test
using BiasedRandomKeyGA

const do_basic_test = true
const do_pathrelinking_test = true
const do_makespan_test = true


const MAXITER = 100000

do_basic_test && include("basictest.jl")
do_pathrelinking_test && include("basicpathrelinking.jl")
do_makespan_test && include("makespan.jl")

