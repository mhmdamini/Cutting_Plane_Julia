function unittest()
    model=Model(HiGHS.Optimizer);
    @variable(model, x[1:2] >= 0, Int)  
    @objective(model, Min, x[1] + 2*x[2])
    @constraint(model,C1, -2*x[1] - 2*x[2]<= -0.5)
    @constraint(model,C2, 2*x[1] - x[2] <= 1.5)
    @constraint(model,C3, x[1]  <= 1)
    @constraint(model,C4, x[2]<= 1)
    optimize_model(model,OPT,false)
    print(model)
    return model
end


