using DataFrames

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