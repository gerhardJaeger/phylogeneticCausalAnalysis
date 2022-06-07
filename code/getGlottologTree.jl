cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

##

using ProgressMeter

##



ENV["PYTHON"] = readchomp(`which python`)
Pkg.build("PyCall")
using PyCall

ete3 = pyimport("ete3")
##

glottologF = "../data/tree_glottolog_newick.txt"

isfile(glottologF) || download(
    "https://cdstar.eva.mpg.de/bitstreams/EAEA0-F8BB-0AB6-96FA-0/tree_glottolog_newick.txt",
    glottologF
)

##

raw = readlines(glottologF);

##

trees = []

for ln in raw
    ln = strip(ln)
    ln = replace(ln, r"\'[A-ZÄÖÜ][^[]*\[" => "[")
    ln = replace(ln, r"\][^']*\'" => "]")
    ln = replace(ln, r"\[|\]" => "")
    ln = replace(ln, ":1" => "")
    push!(
        trees,
        ete3.Tree(ln, format=1)
    )
end

##

glot = ete3.Tree()

for t in trees
    glot.add_child(t)
end

nonLeaves = [nd.name for nd in glot.traverse()
             if (nd.name != "") & !nd.is_leaf()
]

@showprogress for nm in nonLeaves
    nd = (glot & nm)
    nd.name = ""
    nd.add_child(name=nm)
end


##

open("../data/glottolog.tre", "w") do io
    write(io, glot.write(format=9))
end

