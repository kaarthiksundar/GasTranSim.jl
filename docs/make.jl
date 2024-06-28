using Documenter, GasTranSim

makedocs(
    sitename = "GasTranSim",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    authors = "Shriram Srinivasan and Kaarthik Sundar",
    pages = [
        "Introduction" => "index.md",
        "Theory" => "theory.md",
        "Input Files" => "inputs.md",
        "Simulation Engine" => "simulation.md",
        "Examples" => "examples.md"
    ],
)

deploydocs(repo = "https://github.com/kaarthiksundar/GasTranSim.jl.git",
    versions = ["stable" => "v^", "v#.#"])
