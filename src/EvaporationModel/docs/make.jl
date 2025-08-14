using Documenter, EvaporationModel

makedocs(;
    sitename="EvaporationModel",
    pages=["Home" => "index.md"],
    modules=[EvaporationModel],
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/olivierbonte/DifferentiableEvaporation.git", push_preview=true
)
