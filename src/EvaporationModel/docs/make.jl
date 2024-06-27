using Documenter, EvaporationModel

makedocs(;
    sitename="EvaporationModel",
    pages = [
        "Home" => "index.md",
        hide("soil.md")
    ]
)