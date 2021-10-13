using Documenter, PathIntegrationMethod

makedocs(
    modules = [PathIntegrationMethod],
    sitename = "PathIntegrationMethod.jl",
    authors = "Henrik T Sykora",
    clean = false,
    # format = format,
    # checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Path Integration Method" => 
        ["PathIntegration" => "pathintegrationmethod.md",
        # "Specific SDE classes" => "sdeclasses.md",
        # "Interpolation" => "interpolation.md",
        # "Time stepping" => "timestepping.md",
        # "Integration" => "integration.md"
        ],
        # "Examples" => "examples.md"
    ]
)

