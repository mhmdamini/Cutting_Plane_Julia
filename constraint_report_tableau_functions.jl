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
function Generate_optimal_simplex_tableau(norm_model)
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

function simplexes(model,optimizer)
    normal_model = create_normal_form_model(model)
    optimize_model(normal_model, optimizer, false)
    first_simplex_tableau = generate_first_simplex_tableau(normal_model)
    optimal_simplex_tableau, reduced_size_tableau, basis_vars, all_vars, structural_vars = generate_optimal_simplex_tableau(normal_model)
    return first_simplex_tableau,optimal_simplex_tableau
end


function generate_optimal_simplex_tableau(norm_model)
    c1, A, constant, b_upper, b_lower, var_bounds, vars, x_lower, x_upper, cons_list, sense = get_lp_data(norm_model)
    structural_vars = vcat(vars[1:size(A, 2) - size(A, 1)], "z")
    vtype, basis_vars = report_basis(norm_model)
    basic_constraint = get_constraint_status(norm_model)
    if !isempty(basic_constraint)
        basic_constraint_order = map(x -> x isa String ? parse(Int, x) : x, basic_constraint)
        slacks_in_basis = [vars[length(structural_vars) - 1 + basic_constraint_order[i]] for i in 1:length(basic_constraint)]
        basis_vars = vcat(basis_vars, slacks_in_basis)
    end
    df = DataFrame(Matrix(A), Symbol.(vars));
    B = Matrix(df[:, Symbol.(basis_vars)]);
    inv_B = inv(B);
    BB = inv_B * A
    tol = 1e-3
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
    missing_vars = setdiff(structural_vars, reduced_size_df.BasisVars)
    if !isempty(missing_vars)    
        num_other_cols = ncol(reduced_size_df) - 1  # Number of other columns except 'BasisVars'
        zeros_matrix = zeros(length(missing_vars), num_other_cols)
        missing_df = DataFrame(BasisVars=missing_vars)
        for i in 1:num_other_cols
            column_name = names(reduced_size_df)[i+1]
            missing_df[!, column_name] = zeros_matrix[:, i]
        end
        combined_df = vcat(reduced_size_df, missing_df)
        combined_df[indexin(structural_vars, combined_df.BasisVars),:]
        return optimal_simplex_tableau, combined_df, basis_vars, vars, structural_vars
    else
        return optimal_simplex_tableau, reduced_size_df, basis_vars, vars, structural_vars
    end
end