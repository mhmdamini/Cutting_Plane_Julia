using DataFrames
using Base: Float64
using JuMP

#This function processes a given vector of floating-point numbers and identifies the values that have a non-zero fractional part. It returns two vectors: one containing the fractional values and the other containing their corresponding indices in the original vector.
"""function find_fractional_elements(vector::Vector{Float64})::Tuple{Vector{Float64}, Vector{Int}}
    # Vectors to store fractional values and their indices
    fractional_values = Float64[]
    fractional_indices = Int[]

    # Iterate through the input vector to identify fractional elements
    for (index, value) in enumerate(vector)
        # Skip integer values
        if isinteger(value)
            continue
        end

        # Check if the fractional part is non-zero
        fractional_part = value - floor(value)
        if fractional_part != 0.0
            push!(fractional_values, value)
            push!(fractional_indices, index)
        end
    end

    # Return the vectors containing fractional values and their indices
    return fractional_values, fractional_indices
end"""
# Function to find fractional elements in a vector of real numbers
function find_fractional_elements(vector::AbstractVector{<:Real})::Tuple{Vector{Float64}, Vector{Int}}
    # Vectors to store fractional values and their indices
    fractional_values = Float64[]
    fractional_indices = Int[]

    # Iterate through the input vector to identify fractional elements
    for (index, value) in enumerate(vector)
        # Convert value to Float64 for calculations if not already
        value_float = Float64(value)

        # Skip integer values
        if isinteger(value_float)
            continue
        end

        # Check if the fractional part is non-zero
        fractional_part = value_float - floor(value_float)
        if fractional_part != 0.0
            push!(fractional_values, value_float)
            push!(fractional_indices, index)
        end
    end

    # Return the vectors containing fractional values and their indices
    return fractional_values, fractional_indices
end

# Function to add a constraint based on the fractional value of a variable
"""function generate_subproblem_with_split(relaxed_model, fractional_var_index, constraint_direction::Int, constraint_name::String)
    # Copy the input model to avoid modifying the original one
    subproblem_model = copy(relaxed_model)

    # Optimize the copied model
    optimize_model(subproblem_model, OPT, false)

    # Get variable values after optimization
    variable_values, _ = get_variable_values(subproblem_model)

    # Find fractional elements
    fractional_value, fractional_indices = find_fractional_elements(variable_values)

    # Get the specific variable to constrain
    fractional_var = all_variables(subproblem_model)[fractional_indices[fractional_var_index]]

    # Calculate the RHS value for the constraints
    rhs_floor = floor(fractional_value[fractional_var_index])

    # Add the appropriate constraint based on the 'constraint_direction' argument
    if constraint_direction == 1
        new_constraint = @constraint(subproblem_model, fractional_var <= rhs_floor)
    else
        new_constraint = @constraint(subproblem_model, -fractional_var <= -(rhs_floor + 1))
    end

    # Set the name of the new constraint for easy identification
    set_name(new_constraint, constraint_name)

    return subproblem_model
end"""


function generate_subproblem_with_split(relaxed_model::Model, fractional_var_index::Int, constraint_direction::Int, constraint_name::String)
    # Copy the input model to avoid modifying the original one
    subproblem_model = copy(relaxed_model)

    # Optimize the copied model
    optimize_model(subproblem_model, OPT, false)

    # Get variable values after optimization
    variable_values_any, _ = get_variable_values(subproblem_model)
    
    # Ensure all values are of type Float64
    variable_values = Float64.(variable_values_any)

    # Find fractional elements
    fractional_values, fractional_indices = find_fractional_elements(variable_values)

    # Check if the fractional_var_index is within bounds
    if fractional_var_index > length(fractional_values)
        error("fractional_var_index out of bounds. No such fractional element.")
    end

    # Get the specific variable to constrain
    fractional_var = all_variables(subproblem_model)[fractional_indices[fractional_var_index]]

    # Calculate the RHS value for the constraints
    rhs_floor = floor(fractional_values[fractional_var_index])

    # Add the appropriate constraint based on the 'constraint_direction' argument
    if constraint_direction == 1
        new_constraint = @constraint(subproblem_model, fractional_var <= rhs_floor)
    else
        new_constraint = @constraint(subproblem_model, -1*fractional_var <= -(rhs_floor + 1))
    end

    # Set the name of the new constraint for easy identification
    set_name(new_constraint, constraint_name)

    return subproblem_model
end

# Function to generate subproblems based on the disjunction for a specified variable index
function create_disjunctive_subproblems(relaxed_model, fractional_var_index, constraint_names::Vector{String})
    if length(constraint_names) != 2
        error("constraint_names must contain exactly two names.")
    end

    # Generate the subproblem with the first disjunction (<= floor)
    subproblem1 = generate_subproblem_with_split(relaxed_model, fractional_var_index, 1, constraint_names[1])

    # Generate the subproblem with the second disjunction (>= ceil)
    subproblem2 = generate_subproblem_with_split(relaxed_model, fractional_var_index, 2, constraint_names[2])

    return subproblem1, subproblem2
end


# Function to compute the rays of the polyhedron described by the optimal simplex tableau
function compute_rays(df, structural_vars, non_basic_vars, all_vars)
    missing_vars = setdiff(structural_vars, df.BasisVars)  
    zeros_matrix = zeros(length(missing_vars), size(df)[2] - 1)
    
    var_names = vcat("BasisVars", all_vars,"rhs")
    missing_df = DataFrame(hcat(missing_vars, zeros_matrix), Symbol.(var_names))
    combined_df = vcat(df, missing_df)
    
    var_order_dict = Dict(structural_vars .=> 1:length(structural_vars))
    combined_df.Number = getindex.(Ref(var_order_dict), combined_df.BasisVars[1:end])
    sorted_df = sort(combined_df, :Number, rev=false)
    
    rays_df = select(sorted_df, Symbol.(vcat("BasisVars", non_basic_vars)))
    rays_matrix = -Matrix(rays_df[1:end-1, 2:end])
    return rays_matrix, rays_df
end

function solve_linear_problem(ω, points, rays_matrices)
    # Initialize the model with HiGHS optimizer
    model = Model(HiGHS.Optimizer)

    # Number of variables in α
    n = length(ω)

    # Define variables α and β
    @variable(model, α[1:n])
    @variable(model, β==1)

    # Define the objective function: minimize α^T ω
    @objective(model, Min, sum(α[i] * ω[i] for i in 1:n))
    # Define constraints
    @constraint(model, sum(α[i] * points[1][i] for i in 1:n) >= β)
    @constraint(model, sum(α[i] * points[2][i] for i in 1:n) >= β)

    # Constraints for rays_matrices[1]
    for j in 1:size(rays_matrices[1], 2)
        @constraint(model, sum(α[i] * rays_matrices[1][i, j] for i in 1:n) >= 0)
    end

    # Constraints for rays_matrices[2]
    for j in 1:size(rays_matrices[2], 2)
        @constraint(model, sum(α[i] * rays_matrices[2][i, j] for i in 1:n) >= 0)
    end
    optimize_model(model,OPT,false)
    
    # Extract the values of α
    α_values = value.(α)

    # Extract the value of β
    β_value = value(β)
    return model,α_values,β_value
    # Solve the model
    """optimize!(model)

    # Check the solution status
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found.")
        optimal_α = value.(α)
        optimal_β = value(β)
        return optimal_α, optimal_β
    else
        println("The optimization problem did not solve to optimality.")
        return nothing
    end"""
end 