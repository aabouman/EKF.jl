module EKF

using LinearAlgebra
using ForwardDiff: jacobian

export ExtendedKalmanFilter, estimateState

#
struct ExtendedKalmanFilter{T}
    Q::AbstractArray{T}     # Dynamics Noise Covariance
    R::AbstractArray{T}     # Measurement Noise Covariance
    process::Function          # Dynamics Function
    measure::Function       # Measurement Function
    # ∂f_∂x::Function
    # ∂h_∂x::Function

    function ExtendedKalmanFilter(Q::AbstractArray{T}, R::AbstractArray{T},
                                  process::Function, measure::Function) where T
        issymmetric(Q) || throw(ArgumentError("Dynamics noise covariance Matrix, Q, must be symmetric."))
        issymmetric(R) || throw(ArgumentError("Measurement noise covariance Matrix, R, must be symmetric."))

        # ∂f_∂x(est_state, input) = jacobian(state->process(state, input),
        #                                    est_state)
        # ∂h_∂x(est_state) = jacobian(state->ekf.measure(state), est_state)

        new{T}(Q, R, process, measure,
               # ∂f_∂x, ∂h_∂x
               )
    end
end


#
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

    # Fₖ₋₁ = ekf.∂f_∂x(x̂ₖ₋₁⁺, uₖ₋₁)
    Fₖ₋₁ = ∂f_∂x(x̂ₖ₋₁⁺, uₖ₋₁)
    x̂ₖ⁻ = ekf.process(x̂ₖ₋₁⁺, uₖ₋₁)

    Pₖ⁻ = Fₖ₋₁ * Pₖ₋₁⁺ * transpose(Fₖ₋₁) + ekf.Q

    # Hₖ = ekf.∂h_∂x(x̂ₖ⁻)
    Hₖ = ∂h_∂x(x̂ₖ⁻)
    Kₖ = Pₖ⁻ * transpose(Hₖ) * inv(Hₖ * Pₖ⁻ * transpose(Hₖ) + ekf.R)

    ỹₖ = zₖ - ekf.measure(x̂ₖ⁻)
    x̂ₖ⁺ = x̂ₖ⁻ + Kₖ * ỹₖ

    Pₖ⁺ = (I(n) - Kₖ * Hₖ) * Pₖ⁻

    return x̂ₖ⁺, Pₖ⁺
end


function simulate(initState::AbstractArray, initEstimate::AbstractArray,
                  errorCov::AbstractArray, inputs::AbstractArray,
                  numSteps, ekf::ExtendedKalmanFilter)
    state = initState
    est_state = initEstimate

    for t = 1:numSteps
        state = ekf.process(state, inputs[t])
        measurement = ekf.measure(state, inputs[t])
        est_state, errorCov = estimate(est_state, [ext_wrench[t]],
                                       measurement, errorCov, ekf)
        θ̂s[t] = est_state[1];
    end
end

end # module
