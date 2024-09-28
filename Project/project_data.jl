using SparseArrays

nodes_with_generators = [[], [1, 2, 3], [4], [5], [6], [], [7], [], [8, 9], [], []]
nodes_with_consumers = [[1], [], [], [2], [], [3], [], [4], [5], [6], [7]]

num_nodes = 11

generator_nodes = [2, 2, 2, 3, 4, 5, 7, 9, 9]
generator_capacity = [0.02, 0.15, 0.08, 0.07, 0.04, 0.17, 0.17, 0.26, 0.05]
generator_cost = [175, 100, 150, 150, 300, 350, 400, 300, 200]
consumerDemand = [0.10, 0.19, 0.11, 0.09, 0.21, 0.05, 0.04]

reactivepower_ub = 0.03 .* generator_capacity
reactivepower_lb = -0.03 .* generator_capacity

generator_num = length(generator_nodes)

# Voltage (phase) angles (theta_k) constraint.
theta_ub = pi
theta_lb = -pi

# Voltage amplitudes (v_k)constraint.
voltage_ub = 1.02
voltage_lb = 0.98

edges = [(1, 2), (1, 11), (2, 3), (2, 11), (3, 4), (3, 9), (4, 5), (5, 6), (5, 8), (6, 7), (7, 8), (7, 9), (8, 9), (9, 10), (10, 11)]

b_coef = [-20.1, -22.3, -16.8, -17.2, -11.7, -19.4, -10.8, -12.3, -9.2, -13.9, -8.7, -11.3, -14.7, -13.5, -26.7]

g_coef = [4.12, 5.67, 2.41, 2.78, 1.98, 3.23, 1.59, 1.71, 1.26, 1.11, 1.32, 2.01, 2.41, 2.14, 5.06]

b_kl = sparse(zeros(num_nodes, num_nodes))
g_kl = sparse(zeros(num_nodes, num_nodes))


for (i, (k, l)) in enumerate(edges)
    b_kl[k, l] = b_coef[i]
    b_kl[l, k] = b_coef[i]
    g_kl[k, l] = g_coef[i]
    g_kl[l, k] = g_coef[i]
end

function p_kl(v_k, v_l, theta_k, theta_l, k::Int, l::Int)::Float64
    return (v_k^2.0 * g_kl[k, l] - v_k * v_l * g_kl[k, l]*cos(theta_k - theta_l) - v_k * v_l * b_kl[k, l] * sin(theta_k - theta_l))
end

function q_kl(v_k, v_l, theta_k, theta_l, k::Int, l::Int)::Float64
    return (-v_k^2.0 * b_kl[k, l] + v_k * v_l * b_kl[k, l]*cos(theta_k - theta_l) - v_k * v_l * g_kl[k, l] * sin(theta_k - theta_l))
end
