using Pkg
using JuMP

# Load required packages
function setup_environment()
    Pkg.activate("path/to/your/environment")

    # Install required packages if not already installed
    Pkg.add("Dualization")
    Pkg.add("HiGHS")
end

# Import functions from different files
function load_functions()
    include("libraries_and_basic_functions.jl")
    include("constraint_report_tableau_functions.jl")
    include("form_model_functions.jl")
    include("rays_function.jl")
    include("VPC.jl")
    include("dual_functions.jl")
end

OPT = "HiGHS"

# Main execution function to keep flow clean
function main(model_path)
    setup_environment()
    load_functions()


    try
        model = read_model(model_path)
        optimize_model(model, OPT, false)

        relaxed_model = relax_model(model)
        optimize_model(relaxed_model, OPT, false)

        standard_model = create_standard_form_model(relaxed_model)
        optimize_model(standard_model, OPT, false)

        objective_value_primary, termination_status_value = retrieve_objective_values(standard_model)

        subproblems = create_disjunctive_subproblems(standard_model, 1, ["D1", "D2"])
        optimize_model(subproblems[1], OPT, false)
        optimize_model(subproblems[2], OPT, false)

        result1 = process_subproblem(subproblems[1], OPT)
        result2 = process_subproblem(subproblems[2], OPT)

        return (
            objective_value_primary, 
            termination_status_value, 
            result1, 
            result2,
            get_variable_values(subproblems[1]),
            get_variable_values(subproblems[2])
        )
    catch e
        println("An error occurred: ", e)
        return nothing
    end
end

# Function to process a single subproblem
function process_subproblem(subproblem::Model, optimizer::String)
    norm_model = create_normal_form_model(subproblem)
    optimize_model(norm_model, optimizer, false)

    first_simplex_tableau = generate_first_simplex_tableau(norm_model)
    optimal_simplex_tableau, reduced_size_tableau, basis_vars, all_vars, structural_vars = generate_optimal_simplex_tableau(norm_model)
    non_basic_vars = setdiff(all_vars, basis_vars)
    rays_matrix, rays_df = compute_rays(reduced_size_tableau, structural_vars, non_basic_vars, all_vars)

    return (
        first_simplex_tableau, 
        optimal_simplex_tableau, 
        rays_matrix,
        get_variable_values(norm_model)
    )
end

# Execute the main function
results = main("p0033_presolved.mps.gz")

# Access results if needed. For example:
if results !== nothing
    objective_value_primary, termination_status_value, result1, result2, variables1, variables2 = results

    # result1 and result2 are tuples returned from process_subproblem
    first_simplex_tableau1, optimal_simplex_tableau1, rays_matrix1, variables_norm1 = result1
    first_simplex_tableau2, optimal_simplex_tableau2, rays_matrix2, variables_norm2 = result2
    
    # Now you have access to all the returned results
    println("Results successfully obtained.")
end