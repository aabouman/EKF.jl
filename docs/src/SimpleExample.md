```@meta
EditURL = "<unknown>/examples/SimpleExample.jl"
```

# SimpleExample

PREAMBLE

PKG_SETUP

## Extened Kalman Filter Example

```@example SimpleExample
using LinearAlgebra: Diagonal, I
using EKF
using Random
using Plots
```

### Defining The Process and Measurement Functions

```@example SimpleExample
function process(xᵢ::AbstractVector, uᵢ::AbstractVector)
    Δt = .1
    A = [1 Δt 0  0;
         0  1 0  0;
         0  0 1 Δt;
         0  0 0  1]
    return A * xᵢ
end

function measure(xᵢ::AbstractVector)
    x, ẋ, y, ẏ = xᵢ
    return [sqrt(x^2 + y^2); atan(y, x)]
end
```

### Defining The Covariances

```@example SimpleExample
Q = Matrix(Diagonal([0, .1, 0, .1]))
R = Matrix(Diagonal([50^2, 0.005^2]));
nothing #hide
```

```@example SimpleExample
ekf = ExtendedKalmanFilter(Q, R, process, measure)
```

```@example SimpleExample
Random.seed!(1)
initState = rand(4)
initEstimate = rand(4)
errorCov = I(4) * .001
inputs = zeros(1000, 2)

θ̄s, θ̂s = simulate(initState::AbstractArray, initEstimate::AbstractArray,
                  errorCov, inputs::AbstractArray, ekf);
nothing #hide
```

```@example SimpleExample
p1 = plot(θ̂s[:, 1])
plot!(p1, θ̄s[:, 1])
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

