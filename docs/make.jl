using TrypColonies
using Documenter

DocMeta.setdocmeta!(TrypColonies, :DocTestSetup, :(using TrypColonies); recursive=true)

makedocs(;
    modules=[TrypColonies],
    authors="AndreasKuhn-ak <andreaskuhn92@gmx.net> and contributors",
    sitename="TrypColonies.jl",
    format=Documenter.HTML(;
        canonical="https://AndreasKuhn-ak.github.io/TrypColonies.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AndreasKuhn-ak/TrypColonies.jl",
    devbranch="master",
)
