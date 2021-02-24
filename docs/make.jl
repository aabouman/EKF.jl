using Documenter, EKF


Documenter.makedocs(
    sitename = "EKF.jl",
    source  = "src",
    build   = "build",
    clean   = true,
    doctest = true,
    # repo = "https://github.com/aabouman/EKF.jl",
    modules = [EKF],
    pages = ["Home" => "index.md",
             "Documentation" => "documentation.md"
             ]
)

deploydocs(;
    repo = "github.com/aabouman/EKF.jl",
    branch = "gh-pages",
    devurl = "docs",
)
