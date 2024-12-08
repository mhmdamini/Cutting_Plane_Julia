integer_model = unittest()
print(integer_model)
optimize_model(integer_model, optimizer_name, false);
integer_objective_value, integer_termination_status = retrieve_objective_values(integer_model);
relaxed_model = relax_model(integer_model);
optimize_model(relaxed_model, optimizer_name, false);
integer_vars=lp_matrix_data(integer_model).integers;
is_integer_solution = all(f -> isinteger(value(f)), all_variables(relaxed_model)[integer_vars]);
relaxed_model_value, relaxed_model_variables = get_variable_values(relaxed_model);
if is_integer_solution
    println("Optimal integer solution found:");
    for (var, val) in zip(relaxed_model_variables, relaxed_model_value)
        println("$var = $val");
    end
    primary_objective_value, primary_termination_status = retrieve_objective_values(relaxed_model);
    println("The optimal solution is:",primary_objective_value);
else
    println("Integer solution not found, relaxed model result:");
    for (var, val) in zip(relaxed_model_variables, relaxed_model_value)
        println("$var = $val");
    end
    primary_objective_value, primary_termination_status = retrieve_objective_values(relaxed_model);
    println("The optimal solution is:",primary_objective_value);
    standard_form_model = create_standard_form_model1(relaxed_model);
    optimize_model(standard_form_model, optimizer_name, false);
    primary_objective_value, primary_termination_status = retrieve_objective_values(standard_form_model);
    selected_var_index=parse(Int, readline());
    println("applying V-polyhedral disjunctive cut");
    subproblems,D0 = create_disjunctive_subproblems(standard_form_model, selected_var_index, integer_vars,["Disjunction1", "Disjunction2"]);
    optimize_model(subproblems[1], optimizer_name, false);
    optimize_model(subproblems[2], optimizer_name, false);
    first_tableau1, optimal_tableau1, rays_matrix1, variable_norms1,non_basic1,normal_model1 = process_subproblem(subproblems[1], optimizer_name);
    first_tableau2, optimal_tableau2, rays_matrix2, variable_norms2,non_basic1,normal_model2 = process_subproblem(subproblems[2], optimizer_name);
    points1, _ = get_variable_values(normal_model1);
    points2, _ = get_variable_values(normal_model2);
    ray_matrices = [rays_matrix1, rays_matrix2]
    points = [points1, points2]
    ω, _ = get_variable_values(relaxed_model)
    PRLP,α_values,β_value=solve_prlp_problem(ω, points, ray_matrices)
    print(PRLP)
    added_vpc_cut=add_Cut(relaxed_model,α_values, β_value,true)
    optimize_model(added_vpc_cut,optimizer_name,false)
    print(added_vpc_cut)
    primary_objective_value_before_cut, primary_termination_status = retrieve_objective_values(relaxed_model)
    primary_objective_value_after_cut, primary_termination_status = retrieve_objective_values(added_vpc_cut)
    integrality_gap_value = integrality_gap(primary_objective_value_after_cut, integer_objective_value, primary_objective_value_before_cut)
    println("Optimal solution of the main problem is: ", integer_objective_value)
    println("Optimal solution of the relaxed problem is: ", primary_objective_value_before_cut)
    println("Optimal solution of the problem with cut is: ", primary_objective_value_after_cut)
     #Print the integrality gap and final solution
    println("The integrality gap is: ", integrality_gap_value)
    println("Optimal solution is:")
    print_variable_values(get_variable_values(added_vpc_cut))
    vpc_model_value, vpc_model_variables = get_variable_values(added_vpc_cut);
    is_integer_VPC_solution = all(f -> isinteger(value(f)), all_variables(added_vpc_cut)[integer_vars]);
    if is_integer_VPC_solution
        println("Optimal integer solution found:");
        for (var, val) in zip(vpc_model_variables, vpc_model_value)
            println("$var = $val");
        end
    end
    primary_objective_value, primary_termination_status = retrieve_objective_values(relaxed_model);
    println("The optimal solution is:",primary_objective_value_after_cut);
    draw_box_and_print_cut("Monoidal Strengthened Cut")
    
    U,vars2,cons_list2=preparation(subproblems,α_values,ray_matrices)
    α_tilda= monoidal_str(subproblems,selected_var_index,-D0,U,length(vars2),α_values)
    if α_tilda==α_values
        println("There is no Str Cut because α_tilda = α=", α_values)
    else
        println("The obtained V-polyhedral str cut for ",all_variables(relaxed_model)[selected_var_index]," is :")
        added_svpc_cut=add_Cut(relaxed_model,α_tilda, β_value,true)
        optimize_model(added_svpc_cut,optimizer_name,false)
        svpc_model_value, svpc_model_variables = get_variable_values(added_vpc_cut);
        primary_objective_value_after_SVPCcut, primary_termination_status = retrieve_objective_values(added_svpc_cut)
        integrality_gap_value_SVPC = integrality_gap(primary_objective_value_after_SVPCcut, integer_objective_value, primary_objective_value_before_cut)
        println("The integrality gap is: ", integrality_gap_value)
        println("Optimal solution with SVPC:");
        
        for (var, val) in zip(svpc_model_variables, svpc_model_value)
            println("$var = $val");
        end
    end
end