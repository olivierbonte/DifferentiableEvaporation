using Documenter, EvaporationModel

makedocs(; sitename="EvaporationModel", pages=["Home" => "index.md"])

deploydocs(;
    repo="github.com/olivierbonte/DifferentiableEvaporation.git",
    devbranch="feature/online_docs",
)