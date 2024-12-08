# Import necessary libraries
using Pkg
using JuMP
using MathOptInterface
using LinearAlgebra

const MOI = MathOptInterface

# Activate the current directory's environment to ensure correct package dependencies
Pkg.activate(@__DIR__)

# Add the HiGHS optimization package
Pkg.add("HiGHS")

# Install all dependencies listed in the Project.toml file
Pkg.instantiate()

# Function to read a model from a given file path
function read_model(path)
    model = read_from_file(path)
    return model
end

# Function to relax the integrality constraints in the given model
function relax_model(model)
    relaxed_model = copy(model)
    relax_integrality(relaxed_model)
    return relaxed_model
end

# Function to set the optimizer for the given model and perform optimization
function optimize_model(model, optimizer_name, output=false)
    if optimizer_name == "HiGHS"
        set_optimizer(model, optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => output))
    elseif optimizer_name == "GLPK"
        set_optimizer(model, optimizer_with_attributes(GLPK.Optimizer))
    elseif optimizer_name == "Gurobi"
        set_optimizer(model, optimizer_with_attributes(Gurobi.Optimizer))
    elseif optimizer_name == "Clp"
        set_optimizer(model, Clp.Optimizer)
    else
        println("Unknown optimizer: $optimizer_name")
    end
    optimize!(model)
end

# Function to retrieve variable values if the model converged to an optimal solution
function get_variable_values(model)
    if termination_status(model) == MOI.OPTIMAL
        values = []
        names = []
        for var in all_variables(model)
            push!(values, value(var))
            push!(names, name(var))
        end
        return values, names
    else
        println("Optimization did not converge to an optimal solution.")
        return nothing
    end
end

# Function to extract linear programming data from the given model
function get_lp_data(model)
    data = lp_matrix_data(model)
    return (
        data.c,
        data.A,
        data.c_offset,
        data.b_upper,
        data.b_lower,
        data.variable_constraints,
        data.variables,
        data.x_lower,
        data.x_upper,
        data.affine_constraints,
        data.sense
    )
end

# Function to calculate the integrality gap (cut, integer, relaxed)
function integrality_gap(modelc, modelI, modelr)
    """
    Calculate the integrality gap for three given models.

    Args:
        modelc: The continuous relaxation model after generating the cut.
        modelI: The integer programming model.
        modelr: The reference model (usually the linear relaxation).

    Returns:
        The integrality gap as a percentage.
    """
    

    
    # Calculate the integrality gap as a percentage
    gap = (modelc - modelr) / (modelI - modelr) * 100
    
    # Return the calculated integrality gap
    return gap
end


"""This is a alternative when order is not important. (No need now)
function integrality_gap(model1, model2, model3)
    # Helper function to determine if the model is relaxed version
    function is_relaxed(model)
        has_integrality_constraint = false
        for v in all_variables(model)
            if is_integer(v)
                has_integrality_constraint = true
                break
            end
        end
        return !has_integrality_constraint
    end

    # Determine the relaxed and continuous models based on constraints
    function determine_models(modela, modelb, modelc)
        num_constraints_a = length(all_constraints(modela,include_variable_in_set_constraints=false))
        num_constraints_b = length(all_constraints(modelb,include_variable_in_set_constraints=false))
        num_constraints_c = length(all_constraints(modelc,include_variable_in_set_constraints=false))
        
        if num_constraints_a == num_constraints_b + 1
            return modelb, modela  # modelb is relaxed, modela has an extra constraint
        elseif num_constraints_b == num_constraints_a + 1
            return modela, modelb  # modela is relaxed, modelb has an extra constraint
        elseif num_constraints_c == num_constraints_a + 1
            return modela, modelc  # modela is relaxed, modelc has an extra constraint
        elseif num_constraints_a == num_constraints_c + 1
            return modelc, modela  # modelc is relaxed, modela has an extra constraint
        elseif num_constraints_b == num_constraints_c + 1
            return modelc, modelb  # modelc is relaxed, modelb has an extra constraint
        elseif num_constraints_c == num_constraints_b + 1
            return modelb, modelc  # modelb is relaxed, modelc has an extra constraint
        else
            throw(ArgumentError("Failed to distinguish between models."))
        end
    end

    if is_relaxed(model1)
        modelr = model1
        modelc, modelI = determine_models(model2, model3, model1)
    elseif is_relaxed(model2)
        modelr = model2
        modelc, modelI = determine_models(model1, model3, model2)
    elseif is_relaxed(model3)
        modelr = model3
        modelc, modelI = determine_models(model1, model2, model3)
    else
        throw(ArgumentError("All inputs seem to be integer models or the identification heuristic failed."))
    end

    # Calculate objective values for each model
    obj_continuous = objective_value(modelc)  # Objective value of the continuous model
    obj_integer = objective_value(modelI)     # Objective value of the integer model
    obj_relaxed = objective_value(modelr)     # Objective value of the relaxed model

    # Calculate the integrality gap using the formula
    gap = (obj_continuous - obj_relaxed) / (obj_integer - obj_relaxed) * 100

    # Return the computed integrality gap percentage
    return gap
end"""

# Function to retrieve objective values and termination status of an optimization model
function retrieve_objective_values(model::Model)
    # Obtain the primary objective value of the model.
    objective_value_primary = objective_value(model)
    
    # Obtain the termination status of the optimization process.
    termination_status_value = termination_status(model)
    
    # Return the primary objective value and termination status.
    return objective_value_primary, termination_status_value
end

# Function to Print Variable Values
function print_variable_values(variable_values)
    if variable_values !== nothing
        values, names = variable_values
        println("Variable values:")
        for (name, value) in zip(names, values)
            println("  $name = $value")
        end
    else
        println("No variable values to display.")
    end
end

#Function to Extract List of Variables
function list_of_all_variables(model)
    all_vars=all_variables(model)
    return all_vars
end

function draw_box_and_print_cut(message)
    # Define the width of the box based on the length of the message
    box_width = length(message) + 4

    # Draw the box
    println("╔" * "═"^box_width * "╗")
    println("║ " * message * "     ")
    println("╚" * "═"^box_width * "╝")
end