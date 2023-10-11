
using Documenter, Literate, TODOTODO

examples =
    sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

for example in examples
    Literate.markdown(
        example,
        joinpath(@__DIR__, "src"),
        repo_root_url = "https://github.com/TODO/TODOTODO.jl/blob/master",
    )
end

example_mds = first.(splitext.(basename.(examples))) .* ".md"

makedocs(
    modules = [ConstraintTrees],
    clean = false,
    format = Documenter.HTML(
        ansicolor = true,
        canonical = "https://todo.github.io/TODOTODO.jl/stable/",
    ),
    sitename = "TODOTODO.jl",
    linkcheck = false,
    pages = ["README" => "index.md"; example_mds; "Reference" => "reference.md"],
)

deploydocs(
    repo = "github.com/TODO/TODOTODO.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
