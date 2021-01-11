using LinearAlgebra
using ForwardDiff: jacobian


Iz1 = 1/10      # Inertia of the first link
l₁ = 1.0        # Length of the first link
r₁ = 1/2 * l₁   # Length from joint to COM
m₁ = 2.0        # Mass of first link


function InertiaMatrix(θ::AbstractVector{T}) where T
    retmat = [Iz1 + m₁ * r₁^2];
    return retmat
end


function CoriolisMatrix(θ::AbstractVector{T}, θ̇::AbstractVector{T}) where T
    ∇M = jacobian(InertiaMatrix, θ);
    ∇M = reshape(∇M, (length(θ), length(θ), length(θ)));
    C = zeros(T, length(θ), length(θ));

    for i in 1:length(θ), j in 1:length(θ)
        C[i, j] = sum([1/2 * (∇M[k,i,j] + ∇M[j,i,k] - ∇M[i,k,j]) * θ̇[k]
                       for k in length(θ)])
    end
    return C
end


function forwardDynamics(state::AbstractVector, ext_wrench::AbstractVector)
    n = length(state) ÷ 2
    θ = state[1:n]; θ̇ = state[n+1:end]

    M_mat = InertiaMatrix(θ)
    C_mat = CoriolisMatrix(θ, θ̇)
    K_mat = [0.];
    D_mat = [0.];

    mat1 = [zeros(n, n)     Diagonal(ones(n));
            -M_mat \ K_mat  -M_mat \ (C_mat + D_mat)];
    mat2 = [zeros(n, n); 1/M_mat];

    state_dot = mat1 * state + mat2 * ext_wrench;

    return state_dot
end


# function measurementFunction(state::AbstractVector{T},
#                              ext_wrench::AbstractVector) where T
#     n = length(state) ÷ 2
#     θ = state[1:n]; θ̇ = state[n+1:end]
#     return θ
# end


function measurementFunction(state::AbstractVector{T}) where T
    n = length(state) ÷ 2
    θ = state[1:n]; θ̇ = state[n+1:end]

    return [θ[1]]
end
