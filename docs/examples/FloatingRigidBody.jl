import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate();

# %%
# Instantiate correct pacakge
using Revise
using LinearAlgebra: Diagonal, I
using EKF
using StaticArrays
using Rotations
using Plots

# %%
### Defining The Process and Measurement Functions
mass = 1.
inertia = [1. 0 0; 0 1 0; 0 0 1]

"""
Convert a 3-Vector into a skew symmetric matrix.
"""
function vec2skew(ω)
    return ([0    -ω[3]  ω[2];
             ω[3]  0    -ω[1];
            -ω[2]  ω[1]  0]);
end

# The Process Model
function process(x⃗ₖ::Vector, u⃗ₖ::Vector; δt::Real=.01)::Vector
    # Extract components
    p⃗ₖ = x⃗ₖ[1:3]; q⃗ₖ = UnitQuaternion(x⃗ₖ[4:7]);
    v⃗ₖ = x⃗ₖ[8:10]; ω⃗ₖ = x⃗ₖ[11:13]
    f⃗ₖ = u⃗ₖ[1:3]; τ⃗ₖ = u⃗ₖ[4:6]

    # Integrate position
    p⃗ₖ₊₁ = p⃗ₖ + δt * v⃗ₖ
    # Integrate orientation Kinematic step of quaternion, UnitQuaternion
    # constructor renormalizes, removing euler error accumulation
    q⃗ₖ₊₁ = Rotations.params(UnitQuaternion((Rotations.params(q⃗ₖ) + δt * Rotations.kinematics(q⃗ₖ, ω⃗ₖ))...))

    # Build useful mass/coriolis matricies
    M = [mass*I  zeros(3,3); zeros(3,3)  inertia];
    C = [mass*vec2skew(ω⃗ₖ)  zeros(3,3);
         zeros(3,3)  vec2skew(ω⃗ₖ)*inertia];

    temp = [v⃗ₖ; ω⃗ₖ] + δt * inv(M) * ([f⃗ₖ; τ⃗ₖ] - C * [v⃗ₖ; ω⃗ₖ])
    v⃗ₖ₊₁ = temp[1:3]
    ω⃗ₖ₊₁ = temp[4:6]

    x⃗ₖ₊₁ = vcat(p⃗ₖ₊₁, q⃗ₖ₊₁, v⃗ₖ₊₁, ω⃗ₖ₊₁)
    return x⃗ₖ₊₁
end


# The Measurement Model
function measure(x⃗ₖ::Vector)::Vector
    p⃗ₖ = x⃗ₖ[1:3]; q⃗ₖ = UnitQuaternion(x⃗ₖ[4:7]);
    v⃗ₖ = x⃗ₖ[8:10]; ω⃗ₖ = x⃗ₖ[11:13]

    q⃗ₖ = Rotations.params(q⃗ₖ)
    z⃗ₖ = vcat(p⃗ₖ, q⃗ₖ)
    return z⃗ₖ
end

# %%
# ### Defining The Covariances
Q = Matrix(Diagonal([.01 for _ in 1:13]))
R = Matrix(Diagonal([.05 for _ in 1:7]));

# %%
ekf = ExtendedKalmanFilter(Q, R, process, measure);

# %%
initState = [0., 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
initEstimate = rand(13)
errorCov = Matrix(I(13) * .001)
inputs = zeros(1000, 7)

x⃗s, x̂s = simulate(initState, initEstimate, errorCov, inputs, ekf);

# %%
x⃗s


# %%
p1 = plot(x̂s[:, 1])
plot!(p1, x⃗s[:, 1])
