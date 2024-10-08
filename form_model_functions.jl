using JuMP
using MathOptInterface

# Function to check the optimization sense (minimization or maximization) of a given model

function get_objective_sense_description(model::Model)::String
    # Obtain the objective sense (minimization or maximization) from the model
    objective_sense = objective_sense(model)
    
    # Translate the objective sense to a human-readable string
    if objective_sense == MathOptInterface.MIN_SENSE
        return "Minimization"
    elseif objective_sense == MathOptInterface.MAX_SENSE
        return "Maximization"
    else
        return "Unknown"
    end
end

# Function to create a standard form model from the given model
"""function create_standard_form_model(model)
    c, A, constant, b_upper, b_lower, var_bounds, vars, x_lower, x_upper, cons_list, sense = get_lp_data(model)
    std_model = Model(HiGHS.Optimizer)
    n = length(all_variables(model))
    
    # Define the decision variables
    @variable(std_model, α[1:n] >= 0)
    @constraint(std_model, A * α .<= b_upper)
    
    # Define lower bound constraints
    for i in 1:length(x_lower)
        if x_lower[i] != 0
            @constraint(std_model, -α[i] <= -x_lower[i])
        end
    end
    
    # Define upper bound constraints
    @constraint(std_model, I(n) * α .<= x_upper)
    
    # Set the objective function
    @objective(std_model, Max, -sum(c[i] * α[i] for i in 1:n) - constant)
    
    return std_model
end"""

# Function to create the normal form model from the given model
function create_normal_form_model(model)
    c1, A, constant, b_upper, b_lower, var_bounds, vars, x_lower, x_upper, cons_list, sense = get_lp_data(model)
    B = Matrix(A)
    n = length(all_variables(model))
    I_matrix = I(length(cons_list))
    BI = hcat(B, I_matrix)
    ns = size(BI)[1]
    
    c = vcat(c1, zeros(ns - n))
    norm_model = Model(HiGHS.Optimizer)
    @variable(norm_model, β[1:n] >= 0)
    @variable(norm_model, s[1:ns] >= 0)
    
    concatenated_vars = vcat(β, s)
    @constraint(norm_model, BI * concatenated_vars .== b_upper)
    @objective(norm_model, Max, sum(c[i] * concatenated_vars[i] for i in 1:ns) + constant)
    
    return norm_model
end


function create_standard_form_model(model)
    c, A, constant, b_upper, b_lower, var_bounds, vars, x_lower, x_upper, cons_list, sense = get_lp_data(model)
    std_model = Model(HiGHS.Optimizer)
    n = length(all_variables(model))
    
    # Define the decision variables; no need for >= 0 since we'll handle bounds in constraints
    @variable(std_model, α[1:n])
    
    # Inequality constraints for Ax <= b_upper
    @constraint(std_model, A * α .<= b_upper)
    
    # Define lower bound constraints
    for i in 1:length(x_lower)
        if x_lower[i] != -Inf # Handle -Inf by not adding a specific constraint
            @constraint(std_model, -α[i] <= -x_lower[i])
        end
    end

    # Define upper bound constraints
    for i in 1:length(x_upper)
        if x_upper[i] != Inf # Handle Inf by not adding a specific constraint
            @constraint(std_model, α[i] <= x_upper[i])
        end
    end
    
    # Set the objective function (including constant term properly)
    @objective(std_model, Max, -sum(c[i] * α[i] for i in 1:n) - constant)
    
    return std_model
end
