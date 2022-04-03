
using JuMP, Juniper, Ipopt, DataFrames
using CSV
using Dates
import XLSX
using Interpolations
using Plots
using LinearAlgebra

pyplot() #utilizando o pyplot

#t1 = now()
optimizer = Juniper.Optimizer
nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
model = Model(optimizer_with_attributes(optimizer, "nl_solver" => nl_solver))#,"mip_solver"=>mip_solver))

usinas = 3

@variable(model, x[i = 1:usinas] ≥ 0)

@objective(model, Min, sum(x[i] for i = 1:usinas))#REVER tem que multiplicar o deficit por alguma coisa!!

for k = 1:2
    @constraint(model,[i = k],
            x[i]^k - 1 == 0)
end

optimize!(model)


if objective_value(model) > 0
    println("Programa rodou")
end
