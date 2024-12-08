using Dualization
using JuMP

# Function to convert digits to their Unicode subscript representations
function to_unicode_subscript(n::Int, subscripts::Vector{String})::String
    digits_str = string(n)
    return join(subscripts[parse(Int, digit) + 1] for digit in digits_str)
end

# Generate a list of Unicode subscript representations from 1 to n with a given prefix
function generate_unicode_subscript_list(n::Int, prefix::String)::Vector{String}
    subscripts = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"]
    return [prefix * to_unicode_subscript(i, subscripts) for i in 1:n]
end

# Dualize a given JuMP model and rename its variables with Unicode subscripts
function create_dual_model_with_subscripts(model::Model)::Model
    dual_model = dualize(model)
    num_variables = length(all_variables(dual_model))

    # Generate subscript names for dual model variables
    variable_subscript_names = generate_unicode_subscript_list(num_variables, "d")

    # Rename each variable in the dual model with the corresponding subscript name
    for (index, variable) in enumerate(all_variables(dual_model))
        set_name(variable, variable_subscript_names[index])
    end

    # Print the number of variables in the dual model
    println("Number of variables in dual model: ", num_variables)

    return dual_model
end

# Example usage: Define and dualize a simple JuMP model
#Test
#model = Model()
#@variable(model, x >= 0)
#@constraint(model, con, x >= 1)
#@objective(model, Min, x)
#dual_model = create_dual_model_with_subscripts(model)