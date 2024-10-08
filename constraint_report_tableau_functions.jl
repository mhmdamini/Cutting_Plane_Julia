using MathOptInterface

const MOI = MathOptInterface

# Function to get the status of constraints in the optimal solution
function get_constraint_status(model)
    status_list = []
    for constraint in all_constraints(model, include_variable_in_set_constraints=false)
        if MOI.get(model, MOI.ConstraintBasisStatus(), constraint) == MOI.BASIC
            push!(status_list, MOI.get(model, MOI.ConstraintName(), constraint))
        end
    end
    return status_list
end

# Function to report the basis status of variables in the model
function report_basis(model)
    variable_basis_status = []
    basis_variables = []
    for i in 1:length(all_variables(model))
        var = all_variables(model)[i]
        status = MOI.get(model, MOI.VariableBasisStatus(), var)
        push!(variable_basis_status, (var, status))
        if status == MOI.BASIC
            push!(basis_variables, var)
        end
    end
    return variable_basis_status, basis_variables
end

using DataFrames

# Function to generate the initial simplex tableau for a given normalized model
function generate_first_simplex_tableau(norm_model)
    c1, A, constant, b_upper, b_lower, var_bounds, vars, x_lower, x_upper, cons_list, sense = get_lp_data(norm_model)
    
    slack_vars = vcat(vars[size(A, 2) - size(A, 1) + 1:end], "z")
    rhs = vcat(b_upper, 0)
    var_names = vcat(vars, "rhs")
    
    initial_simplex_tableau = DataFrame(hcat(vcat(Matrix(A), -c1'), rhs), Symbol.(var_names))
    initial_simplex_tableau[!,:BasisVars] = slack_vars
    select!(initial_simplex_tableau, :BasisVars, names(initial_simplex_tableau)[1:end-1])
    
    return initial_simplex_tableau
end

# Function to generate the optimal simplex tableau for a given normalized model
function generate_optimal_simplex_tableau(norm_model)
    c1, A, constant, b_upper, b_lower, var_bounds, vars, x_lower, x_upper, cons_list, sense = get_lp_data(norm_model)
    
    structural_vars = vcat(vars[1:size(A, 2) - size(A, 1)], "z")
    vtype, basis_vars = report_basis(norm_model)
    slacks_in_basis = get_constraint_status(norm_model)
    
    basis_vars = vcat(basis_vars, slacks_in_basis)
    df = DataFrame(Matrix(A), Symbol.(vars))
    B = Matrix(df[:, Symbol.(basis_vars)])
    inv_B = inv(B)
    BB = inv_B * A
    tol = 1e-3
    #BB[BB .< tol] .= 0.0
    
    index_of_basis = []
    for col in 1:size(BB, 2)
        column_data = BB[:, col]
        if count(x -> abs(x - 1) < tol, column_data) == 1 && count(x -> abs(x) >= tol && abs(x - 1) >= tol, column_data) == 0
            push!(index_of_basis, col)
        end
    end
    
    c_b = c1[index_of_basis]
    all_indices = 1:length(vars)
    non_basic_indices = setdiff(all_indices, index_of_basis)
    
    reduced_costs = []
    for j in non_basic_indices
        A_j = A[:, j]
        reduced_cost = c1[j] - c_b' * inv_B * A_j
        push!(reduced_costs, reduced_cost)
    end
    
    z_row = zeros(length(c1))
    z_row[non_basic_indices] = reduced_costs
    rhs = inv_B * b_upper
    optimal_value = -c_b' * inv_B * b_lower - constant
    tableau_rhs = vcat(rhs, optimal_value)
    tableau_basic_vars = vcat(basis_vars, "z")
    
    var_names = vcat(vars, "rhs")
    optimal_simplex_tableau = DataFrame(hcat(vcat(BB, z_row'), tableau_rhs), Symbol.(var_names))
    optimal_simplex_tableau[!,:BasisVars] = tableau_basic_vars
    select!(optimal_simplex_tableau, :BasisVars, names(optimal_simplex_tableau)[1:end-1])
    
    reduced_size_df = filter(row -> row[1] in structural_vars, optimal_simplex_tableau)
    return optimal_simplex_tableau, reduced_size_df, basis_vars, vars, structural_vars
end