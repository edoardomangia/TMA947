using JuMP
import Ipopt
using SparseArrays

# Data definitions
generators_at_nodes = [[], [1, 2, 3], [4], [5], [6], [], [7], [], [8, 9], [], []]
consumers_at_nodes = [[1], [], [], [2], [], [3], [], [4], [5], [6], [7]]

total_nodes = 11

# Generator data
generator_nodes = [2, 2, 2, 3, 4, 5, 7, 9, 9]
generator_capacities = [0.02, 0.15, 0.08, 0.07, 0.04, 0.17, 0.17, 0.26, 0.05]
generator_costs = [175, 100, 150, 150, 300, 350, 400, 300, 200]
consumer_demands = [0.10, 0.19, 0.11, 0.09, 0.21, 0.05, 0.04]

# Reactive power limits based on generator capacities
reactive_power_upper_bound = 0.03 .* generator_capacities
reactive_power_lower_bound = -0.03 .* generator_capacities

num_generators = length(generator_nodes)

# Voltage angle constraints
voltage_angle_upper_bound = pi
voltage_angle_lower_bound = -pi

# Voltage magnitude constraints
voltage_magnitude_upper_bound = 1.02
voltage_magnitude_lower_bound = 0.98

# Graph edges representing connections between nodes
edges = [(1, 2), (1, 11), (2, 3), (2, 11), (3, 4), (3, 9), (4, 5), (5, 6), (5, 8), (6, 7), (7, 8), (7, 9), (8, 9), (9, 10), (10, 11)]

# Susceptance and conductance coefficients
susceptance_coefficients = [-20.1, -22.3, -16.8, -17.2, -11.7, -19.4, -10.8, -12.3, -9.2, -13.9, -8.7, -11.3, -14.7, -13.5, -26.7]
conductance_coefficients = [4.12, 5.67, 2.41, 2.78, 1.98, 3.23, 1.59, 1.71, 1.26, 1.11, 1.32, 2.01, 2.41, 2.14, 5.06]

# Initialize conductance and susceptance matrices
susceptance_matrix = sparse(zeros(total_nodes, total_nodes))
conductance_matrix = sparse(zeros(total_nodes, total_nodes))

for (i, (node_k, node_l)) in enumerate(edges)
    susceptance_matrix[node_k, node_l] = susceptance_coefficients[i]
    susceptance_matrix[node_l, node_k] = susceptance_coefficients[i]
    conductance_matrix[node_k, node_l] = conductance_coefficients[i]
    conductance_matrix[node_l, node_k] = conductance_coefficients[i]
end

# Define active and reactive power flow functions
function calculate_active_power(voltage_mag_k, voltage_mag_l, voltage_angle_k, voltage_angle_l, k::Int, l::Int)::Float64
    return (voltage_mag_k^2.0 * conductance_matrix[k, l] - 
            voltage_mag_k * voltage_mag_l * conductance_matrix[k, l] * cos(voltage_angle_k - voltage_angle_l) - 
            voltage_mag_k * voltage_mag_l * susceptance_matrix[k, l] * sin(voltage_angle_k - voltage_angle_l))
end

function calculate_reactive_power(voltage_mag_k, voltage_mag_l, voltage_angle_k, voltage_angle_l, k::Int, l::Int)::Float64
    return (-voltage_mag_k^2.0 * susceptance_matrix[k, l] + 
            voltage_mag_k * voltage_mag_l * susceptance_matrix[k, l] * cos(voltage_angle_k - voltage_angle_l) - 
            voltage_mag_k * voltage_mag_l * conductance_matrix[k, l] * sin(voltage_angle_k - voltage_angle_l))
end

# Create optimization model
model = Model(Ipopt.Optimizer)

# Active power generated by each generator constrained by its capacity
@variable(model, 0 <= active_power[i = 1:num_generators] <= generator_capacities[i])

# Voltage angles for each node, constrained between theta_lb and theta_ub
@variable(model, voltage_angle_lower_bound <= voltage_angle[1:total_nodes] <= voltage_angle_upper_bound)

# Voltage magnitudes for each node, constrained between voltage_lb and voltage_ub
@variable(model, voltage_magnitude_lower_bound <= voltage_mag[1:total_nodes] <= voltage_magnitude_upper_bound)

# Reactive power generated by each generator constrained between reactivepower_lb and reactivepower_ub
@variable(model, reactive_power_lower_bound[i] <= reactive_power[i = 1:num_generators] <= reactive_power_upper_bound[i])

# Minimize total cost: sum of each generator's active power multiplied by its cost
@objective(model, Min, sum(active_power[i] * generator_costs[i] for i in 1:num_generators))

# Constraints
active_power_constraints = []
reactive_power_constraints = []

for node in 1:total_nodes
    # Active power balance constraint
    constraint1 = @NLconstraint(model, 
        - sum(
                voltage_mag[node]^2.0 * conductance_matrix[node, j] - 
                voltage_mag[node] * voltage_mag[j] * conductance_matrix[node, j] * cos(voltage_angle[node] - voltage_angle[j]) - 
                voltage_mag[node] * voltage_mag[j] * susceptance_matrix[node, j] * sin(voltage_angle[node] - voltage_angle[j])
                for j in 1:total_nodes
                )
        + sum(
                active_power[j] for j in generators_at_nodes[node]
                )
        - sum(
                consumer_demands[j] for j in consumers_at_nodes[node]
                )
        == 0)
    push!(active_power_constraints, constraint1)
  
    # Reactive power balance constraint
    constraint2 = @NLconstraint(model,
        - sum(
                -voltage_mag[node]^2.0 * susceptance_matrix[node, j] + 
                voltage_mag[node] * voltage_mag[j] * susceptance_matrix[node, j] * cos(voltage_angle[node] - voltage_angle[j]) - 
                voltage_mag[node] * voltage_mag[j] * conductance_matrix[node, j] * sin(voltage_angle[node] - voltage_angle[j])
                for j in 1:total_nodes
                )
        + sum(
                reactive_power[k] for k in generators_at_nodes[node]
                )
        == 0)
    push!(reactive_power_constraints, constraint2)
end

# Optimize the model
optimize!(model)

# Retrieve solution values
active_power_values = JuMP.value.(active_power)
voltage_magnitude_values = JuMP.value.(voltage_mag)
voltage_angle_values = JuMP.value.(voltage_angle)

# Output results
println("\n--- Optimization Results ---")
println("\nTermination status: ", JuMP.termination_status(model))
println("Objective function value (Total Cost): ", round(JuMP.objective_value(model), digits=6))

# Active power generated
println("\nActive Power Generated by Each Generator:")
for i in 1:num_generators
    println("Generator $i: ", round(active_power_values[i], digits=2), " (Max Capacity: ", round(generator_capacities[i], digits=6), ")")
end

# Reactive power generated
println("\nReactive Power Generated by Each Generator:")
for i in 1:num_generators
    println("Generator $i: ", round(JuMP.value(reactive_power[i]), digits=6), " (Reactive Limits: ±", round(reactive_power_upper_bound[i], digits=6), ")")
end

# Voltage magnitudes and angles
println("\nVoltage Magnitudes and Angles at Each Node:")
for i in 1:total_nodes
    println("Node $i: Voltage Magnitude = ", round(voltage_magnitude_values[i], digits=6), ", Voltage Angle = ", round(voltage_angle_values[i], digits=6))
end

# Dual variables / Lagrange multipliers
println("\nDual Variables / Lagrange Multipliers:")

# Format for Active Power Constraints
println("\nActive Power Constraints:")
for (i, dual_value) in enumerate(round.(JuMP.dual.(active_power_constraints), digits=6))
    println("Generator $i: ", dual_value)
end

# Format for Reactive Power Constraints
println("\nReactive Power Constraints:")
for (i, dual_value) in enumerate(round.(JuMP.dual.(reactive_power_constraints), digits=6))
    println("Generator $i: ", dual_value)
end


# Calculate and display power flows
# println("\n--- Power Flows ---")
active_power_flows = zeros(total_nodes, total_nodes)
reactive_power_flows = zeros(total_nodes, total_nodes)

# Calculate power flows
for (k, l) in edges
    active_power_flows[k, l] = calculate_active_power(voltage_magnitude_values[k], voltage_magnitude_values[l], voltage_angle_values[k], voltage_angle_values[l], k, l)
    reactive_power_flows[k, l] = calculate_reactive_power(voltage_magnitude_values[k], voltage_magnitude_values[l], voltage_angle_values[k], voltage_angle_values[l], k, l)
end

println("\nActive Power Flows:")
for k in 1:total_nodes
    for l in 1:total_nodes
        if active_power_flows[k, l] > 0
            println("Active power flow from node $k to node $l: ", round(active_power_flows[k, l], digits=6))
        end
    end
end

println("\nReactive Power Flows:")
for k in 1:total_nodes
    for l in 1:total_nodes
        if reactive_power_flows[k, l] > 0
            println("Reactive power flow from node $k to node $l: ", round(reactive_power_flows[k, l], digits=6))
        end
    end
end