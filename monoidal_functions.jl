include("libraries_and_basic_functions.jl")

function split_lower_bound(model,ind,optimizer_name)
    copy_model=copy(model)
    vars=all_variables(copy_model);
    @objective(copy_model,Min,vars[ind]);
    optimize_model(copy_model,optimizer_name);
    return objective_value(copy_model)
end

function Δ_finder(D,l)
    if l==D
        Δ=D-(l-1);
    else
        Δ=D-l;
    end
    return Δ
end

function monoidal(alpha,u0,Delta)
    m=zeros(length(alpha))
    while true
        ar = [alpha[i] + u0[i] * Delta[i] * m[i] for i in 1:length(alpha)]
        as = [alpha[i] + u0[i] * Delta[i] * (m[i]+1) for i in 1:length(alpha)]
        mxr = maximum(ar)
        r = argmax(ar)
        mns = minimum(as)
        s = argmin(as)
        if mxr<=mns
            return m
            break
        else
            m[r] -= 1
            m[s] += 1
        end
    end
    return m
end

function monoidal_str(subproblems,selected_var_index,D0,U,num_vars,alpha)
    l = [split_lower_bound(subproblems[t], selected_var_index, optimizer_name) for t in 1:length(subproblems)]
    Δ = [Δ_finder(D0[t],l[t]) for t in 1:length(subproblems)]
    u0 = [last(U[i]) for i in 1:length(U)]
    u = [ [-U[j][i] for j in 1:length(U)] for i in 1:num_vars ]
    monoidal_vector=[monoidal(u[t],u0,Δ) for t in 1:num_vars ]
    str_alpha = zeros(Float64, num_vars)
    for i in 1:num_vars
        str_alpha[i] = maximum([u[i][j] + monoidal_vector[i][j] * u0[j] for j in 1:length(u0)])
    end
    return str_alpha+alpha
end


function actives_constraints(problem)
    first_tableau, optimal_tableau, rays_matrix, variable_norms = process_subproblem(problem, optimizer_name);
    slacks_list=first_tableau.BasisVars
    basis_list=optimal_tableau.BasisVars
    nbasis=setdiff(slacks_list,basis_list)
    active_constrains = findall(x -> x in nbasis, slacks_list)
    return active_constrains
end


function preparation(subproblems,α_values,rays_matrix)
    u_ind=[actives_constraints(subproblems[t]) for t in 1:length(subproblems)]

    s_u=[(α_values')rays_matrix[i] for i in 1:length(rays_matrix)]
    _,_, _, _, _, _, vars2, _, _, cons_list2, _ = get_lp_data(subproblems[1])
    U=[zeros(length(cons_list2)),zeros(length(cons_list2))]
    for j in 1:length(U)
        for i in 1:length(u_ind[j])
            U[j][u_ind[j][i]] = s_u[j][i]
        end
    end
    return U,vars2,cons_list2
end