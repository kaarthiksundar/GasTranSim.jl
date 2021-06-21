using Documenter, GasTranSim

makedocs(
    sitename = "GasTranSim",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        assets=[asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)],
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    strict = true,
    authors = "Shriram Srinivasan and Kaarthik Sundar",
    pages = [
        "Introduction" => "index.md"
    ],
)

deploydocs(repo = "https://github.com/kaarthiksundar/GasTranSim.jl.git")