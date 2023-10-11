
using Documenter, Literate, SBMLFBCModels

examples =
    sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

for example in examples
    Literate.markdown(
        example,
        joinpath(@__DIR__, "src"),
        repo_root_url = "https://github.com/COBREXA/SBMLFBCModels.jl/blob/master",
    )
end

example_mds = first.(splitext.(basename.(examples))) .* ".md"

makedocs(
    modules = [ConstraintTrees],
    clean = false,
    format = Documenter.HTML(
        ansicolor = true,
        canonical = "https://todo.github.io/SBMLFBCModels.jl/stable/",
    ),
    sitename = "SBMLFBCModels.jl",
    linkcheck = false,
    pages = ["README" => "index.md"; example_mds; "Reference" => "reference.md"],
)

deploydocs(
    repo = "github.com/COBREXA/SBMLFBCModels.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
