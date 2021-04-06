module EKF

using LinearAlgebra: inv, I, issymmetric, isposdef
using ForwardDiff: jacobian
using Distributions: MvNormal

export ExtendedKalmanFilter, estimateState, simulate

@doc raw"""
`struct ExtendedKalmanFilter{T}`

Extended Kalman Filter struct. Stores the dynamics and measurement functions as
well as their cooresponding covariances.

# Arguments
- `Q::AbstractArray{T}`: Dynamics noise covariance matrix, must be symmetric
- `R::AbstractArray{T}`: Measurement noise covariance matrix, must be symmetric
- `process::Function`: dynamics function, steps the system forward
- `measure::Function`: measurement function

The `process` and `measure` function have the following forms:

```julia
function process(state::AbstractVector, input::AbstractVector)
    ...
    return new_state
end
```

```julia
function measure(state::AbstractVector)
    ...
    return [θ₁, θ₂]
end
```
"""
struct ExtendedKalmanFilter{T}
    n::Int64                # Number of state variables
    m::Int64                # Number of input variables
    Q::AbstractArray{T}     # Dynamics Noise Covariance
    R::AbstractArray{T}     # Measurement Noise Covariance
    process::Function       # Dynamics Function
    measure::Function       # Measurement Function

    function ExtendedKalmanFilter(Q::AbstractMatrix{T}, R::AbstractMatrix{T},
                                  process::Function, measure::Function) where T
        issymmetric(Q) || throw(ArgumentError("Dynamics noise covariance Matrix, Q, must be symmetric."))
        issymmetric(R) || throw(ArgumentError("Measurement noise covariance Matrix, R, must be symmetric."))

        isposdef(Q) || throw(ArgumentError("Dynamics noise covariance Matrix, Q, must be positive semi-definite."))
        isposdef(R) || throw(ArgumentError("Measurement noise covariance Matrix, R, must be positive semi-definite."))

        n = size(Q)[1]
        m = size(R)[1]

        new{T}(n, m, Q, R, process, measure, )
    end
end


@doc raw"""
`estimateState(est_state::AbstractArray{T}, input::AbstractArray{T},
               measurement::AbstractArray{T}, errorCov::AbstractArray{T},
               ekf::ExtendedKalmanFilter{T}) where T`

Estimate the state of the system specified by `ekf`. Returns the new state at
time step ``i+1``

# Arguments
- `est_state::AbstractArray{T}`: Estimated state at time step ``i``
- `input::AbstractArray{T}`: Control input at time step ``i``
- `measurement::AbstractArray{T}`: Measurment of state at time step ``i``
- `errorCov::AbstractArray{T}`: Covariance error
- `ekf::ExtendedKalmanFilter{T}`: ExtendedKalmanFilter struct specifying the dynamics and process aswell as their covariances

The `process` and `measure` function have the following forms:

```julia
function process(state::AbstractVector, input::AbstractVector)
    ...
    return new_state
end
```

```julia
function measure(state::AbstractVector)
    ...
    return [θ₁, θ₂]
end
```

Both of these functions must be differentiable using the `ForwardDiff` package.
"""
function estimateState(est_state::AbstractArray{T}, input::AbstractArray{T},
                       measurement::AbstractArray{T}, errorCov::AbstractArray{T},
                       ekf::ExtendedKalmanFilter{T}) where T
    ∂f_∂x(est_state, input) = jacobian(state->ekf.process(state, input),
                                       est_state)
    ∂h_∂x(est_state) = jacobian(state->ekf.measure(state), est_state)

    # State and Measurement size
    local n = length(est_state); local m = length(measurement);

    # Relabeling
    x̂ₖ₋₁⁺ = est_state; uₖ₋₁ = input;
    zₖ = measurement; Pₖ₋₁⁺ = errorCov;

    Fₖ₋₁ = ∂f_∂x(x̂ₖ₋₁⁺, uₖ₋₁)
    x̂ₖ⁻ = ekf.process(x̂ₖ₋₁⁺, uₖ₋₁)

    Pₖ⁻ = Fₖ₋₁ * Pₖ₋₁⁺ * transpose(Fₖ₋₁) + ekf.Q

    Hₖ = ∂h_∂x(x̂ₖ⁻)
    Kₖ = Pₖ⁻ * transpose(Hₖ) * inv(Hₖ * Pₖ⁻ * transpose(Hₖ) + ekf.R)

    ỹₖ = zₖ - ekf.measure(x̂ₖ⁻)
    x̂ₖ⁺ = x̂ₖ⁻ + Kₖ * ỹₖ

    Pₖ⁺ = (I(n) - Kₖ * Hₖ) * Pₖ⁻

    return x̂ₖ⁺, Pₖ⁺
end


@doc raw"""
`simulate(initState::AbstractArray, initEstimate::AbstractArray,
          errorCov::AbstractArray, inputs::AbstractArray,
          numSteps::Int64, ekf::ExtendedKalmanFilter)`

Simulates the system specified by `ekf::ExtendedKalmanFilter` over time horizon.

# Arguments
- `initState::AbstractArray`: Inital state at time step `0`
- `initEstimate::AbstractArray`: Inital "guess" of the state at time step `0`
- `errorCov::AbstractArray`: Error covariances at each time step
- `inputs::AbstractArray`: List of inputs for each time step
- `numSteps::Int64`: Covariance error
- `ekf::ExtendedKalmanFilter`: ExtendedKalmanFilter struct specifying the dynamics and process as well as their covariances


"""
function simulate(initState::AbstractArray, initEstimate::AbstractArray,
                  errorCov::AbstractArray, inputs::AbstractArray,
                  ekf::ExtendedKalmanFilter)
    numStates = length(initState)
    numSteps = size(inputs)[1]
    state = initState
    est_state = initEstimate

    θ̄s = zeros(numSteps+1, ekf.n)
    θ̄s[1,:] .= state
    θ̂s = zeros(numSteps, ekf.n)

    # Build disturbance distributions
    process_noise_dist = MvNormal(ekf.Q)
    measure_noise_dist = MvNormal(ekf.R)

    for t = 1:numSteps
        state = ekf.process(state, inputs[t,:]) + reshape(rand(process_noise_dist, 1), ekf.n)
        measurement = ekf.measure(state) + reshape(rand(measure_noise_dist, 1), ekf.m)
        est_state, errorCov = estimateState(est_state, inputs[t,:],
                                            measurement, errorCov, ekf)
        θ̄s[t+1,:] .= state
        θ̂s[t,:] .= est_state[1];
    end

    return (θ̄s, θ̂s)
end

end # module
