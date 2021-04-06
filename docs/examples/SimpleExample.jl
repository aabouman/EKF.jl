# # @__NAME__ %%

# %%
# PREAMBLE

# %%
# PKG_SETUP
import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate();

# %%
# ## Extened Kalman Filter Example
using Revise
using LinearAlgebra: Diagonal, I, isposdef
using EKF
using Random
using Plots

# %%
# ### Defining The Process and Measurement Functions
function process(xᵢ::AbstractVector, uᵢ::AbstractVector)
    Δt = .1
    A = [1  Δt  0   0;
         0   1  0   0;
         0   0  1  Δt;
         0   0  0   1]
    return A * xᵢ
end

function measure(xᵢ::AbstractVector)
    x, ẋ, y, ẏ = xᵢ
    return [sqrt(x^2 + y^2); atan(y, x)]
end

# %%
# ### Defining The Covariances
Q = Matrix(Diagonal([1e-10, .1, 1e-10, .1]))
R = Matrix(Diagonal([.5^2, 0.005^2]));

# %%
ekf = ExtendedKalmanFilter(Q, R, process, measure);

# %%
Random.seed!(1)
initState = rand(4)
initEstimate = rand(4)
errorCov = I(4) * .001
inputs = zeros(1000, 2)

x⃗s, x̂s = simulate(initState, initEstimate, errorCov, inputs, ekf);

# %%
p1 = plot(x̂s[:, 1])
plot!(p1, x⃗s[:, 1])
