# Activate Project Environment and Install Necessary Packages
using Pkg
Pkg.activate("environment")
Pkg.add("Dualization")
Pkg.add("HiGHS")

# Import Libraries
using JuMP
using HiGHS

# Import External Functions
include("libraries_and_basic_functions.jl")
include("constraint_report_tableau_functions.jl")
include("form_model_functions.jl")
include("rays_function.jl")
include("VPC.jl")
include("dual_functions.jl")
include("test.jl")

# Define optimizer name
optimizer_name = "HiGHS"

# Run unit tests and obtain the model
model = unittest()

# Optimize the model with the defined optimizer
optimize_model(model, optimizer_name, false)

# Retrieve the objective values
integer_objective_value, integer_termination_status = retrieve_objective_values(model)

# Relax the model to a continuous version
relaxed_model = relax_model(model)

# Optimize the relaxed model
optimize_model(relaxed_model, optimizer_name, false)

# Create a standard form version of the relaxed model
standard_form_model = create_standard_form_model(relaxed_model)

# Optimize the standard form model
optimize_model(standard_form_model, optimizer_name, false)

# Retrieve objective values for the standard form model
primary_objective_value, primary_termination_status = retrieve_objective_values(standard_form_model)

# Create disjunctive subproblems based on the standard form model
subproblems = create_disjunctive_subproblems(standard_form_model, 1, ["Disjunction1", "Disjunction2"])

# Optimize each subproblem
optimize_model(subproblems[1], optimizer_name, false)
optimize_model(subproblems[2], optimizer_name, false)

# Process each subproblem to retrieve tableau and ray information
first_tableau1, optimal_tableau1, rays_matrix1, variable_norms1 = process_subproblem(subproblems[1], optimizer_name)
first_tableau2, optimal_tableau2, rays_matrix2, variable_norms2 = process_subproblem(subproblems[2], optimizer_name)

# Retrieve variable values from subproblems
points1, _ = get_variable_values(subproblems[1])
points2, _ = get_variable_values(subproblems[2])

# Prepare ray matrices and points for linear problem solving
ray_matrices = [rays_matrix1, rays_matrix2]
points = [points1, points2]

# Get variable values from the relaxed model
ω, _ = get_variable_values(relaxed_model)

# Define an initial ζ vector
ζ = ones(2)

# Solve the linear problem with the given ζ, points, and ray matrices
PRLP, α_values, β_value = solve_linear_problem(ζ, points, ray_matrices)

# Add a new constraint to the relaxed model
@constraint(relaxed_model, α_values' * all_variables(relaxed_model) >= β_value)

# Optimize the relaxed model again with the new constraint
optimize_model(relaxed_model, optimizer_name, false)

# Retrieve the new objective values after adding the cut
primary_objective_value_after_cut, primary_termination_status = retrieve_objective_values(relaxed_model)

# Calculate the integrality gap
integrality_gap_value = integrality_gap(primary_objective_value_after_cut, integer_objective_value, primary_objective_value)


println("Optimal solution of the main problem is: ", integer_objective_value)
println("Optimal solution of the relaxed problem is: ", primary_objective_value)
println("Optimal solution of the problem with cut is: ", primary_objective_value_after_cut)
# Print the integrality gap and final solution
println("The integrality gap is: ", integrality_gap_value)
println("Optimal solution is:")
print_variable_values(get_variable_values(relaxed_model))