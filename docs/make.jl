using Literate, Documenter, EKF

Literate.markdown("examples/SimpleExample.jl", "src/"; documenter=true)

Documenter.makedocs(
    sitename = "EKF.jl",
    source  = "src",
    build   = "build",
    clean   = true,
    doctest = true,
    # repo = "https://github.com/aabouman/EKF.jl",
    modules = [EKF],
    pages = ["Home" => "index.md",
             "Example" => "SimpleExample.md",
             "Documentation" => "documentation.md"
             ]
)

deploydocs(;
    repo = "github.com/aabouman/EKF.jl",
    branch = "gh-pages",
    devurl = "docs",
)
