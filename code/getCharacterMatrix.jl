cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

##
using CSV, DataFrames, ZipFile, Pipe
using ProgressMeter
using DataStructures

##


ENV["PYTHON"] = readchomp(`which python`)
Pkg.build("PyCall")
using PyCall

ete3 = pyimport("ete3")
##

d = CSV.File("../data/data_cues.txt") |> DataFrame
d[19, 1] = "stan1325"


##


asjp19ClusteredF = "../data/asjp19Clustered.csv"

isfile(asjp19ClusteredF) || download(
    "https://osf.io/kxm64/download",
    "../data/asjp19Clustered.csv",
)


asjp19Clustered = CSV.read(asjp19ClusteredF, DataFrame)

##



asjpLanguages = "../data/languages.csv"

asjpZip = "../data/asjp-v19.1.zip"

isfile(asjpLanguages) || begin
    download(
        "https://zenodo.org/record/3843469/files/lexibank/asjp-v19.1.zip?download=1",
        asjpZip,
    )
    r = ZipFile.Reader(asjpZip)
    languagesF = [f for f in r.files if occursin("languages", f.name)] |> first
    open(asjpLanguages, "w") do io
        write(io, read(languagesF, String))
    end
    rm(asjpZip)
end

languages = CSV.File(asjpLanguages) |> DataFrame

##



longnames = [
    replace(join([r.classification_wals, r.Name], "."), "-" => "_") for
    r in eachrow(languages)
]

insertcols!(languages, 1, :longname => longnames)

##

asjp = @pipe languages |>
             select(_, [:longname, :Glottocode]) |>
             dropmissing(_, :Glottocode) |>
             innerjoin(
                 _,
                 asjp19Clustered,
                 on=:longname => :language
             ) |>
             filter(x -> x.Glottocode ∈ d.Glottocode, _)
##


taxa_df = @pipe asjp |>
                unique(_, [:longname, :concept]) |>
                groupby(_, :longname) |>
                combine(_, nrow) |>
                leftjoin(_, asjp, on=:longname) |>
                select(_, [:longname, :Glottocode, :nrow]) |>
                unique |>
                sort(_, :nrow, rev=true) |>
                select(_, [:longname, :Glottocode]) |>
                unique(_, :Glottocode) |>
                innerjoin(
                    _,
                    d[:, [:Glottocode, :Language]],
                    on=:Glottocode,
                )
##

glot = ete3.Tree("../data/glottolog.tre")

##


for nd in glot.get_leaves()
    doculects = filter(x -> x.Glottocode == nd.name, d).Language
    for d in doculects
        nd.add_child(name=d)
    end
end


doculects = intersect(glot.get_leaf_names(), d.Language)

glot.prune(doculects)

##


df = @pipe asjp |>
           select(_, [:Glottocode, :concept, :word]) |>
           groupby(_, [:Glottocode, :concept]) |>
           combine(_, :word => join => :words) |>
           unstack(_, :Glottocode, :concept, :words)

##



concepts = unique(asjp.concept)

sounds = first.(sort(unique(split(join(asjp.word), ""))))

##

scDict = Dict()
for (i, l) in enumerate(df.Glottocode), (j, c) in enumerate(concepts), s in sounds
    if ismissing(df[i, j+1])
        scDict[l, c, s] = "-"
    else
        scDict[l, c, s] = Int(s ∈ df[i, j+1])
    end
end

##


ccc = unique(asjp.cc)

lc = Dict()
for l in taxa_df.Glottocode, c in concepts
    lc[l, c] = false
end

for (l, c) in zip(asjp.Glottocode, asjp.concept)
    lc[l, c] = true
end

cc2c = Dict(zip(asjp.cc, asjp.concept))

ccDict = Dict()
for l in taxa_df.Glottocode, cc in ccc
    if lc[l, cc2c[cc]]
        ccDict[l, cc] = "0"
    else
        ccDict[l, cc] = "-"
    end
end

for (l, cc) in zip(asjp.Glottocode, asjp.cc) |> unique
    ccDict[l, cc] = "1"
end

##


ccMtx = Array{String}(undef, size(taxa_df, 1), size(ccc, 1))
for (i, l) in enumerate(taxa_df.Language), (j, cc) in enumerate(ccc)
    g = taxa_df.Glottocode[i]
    ccMtx[i, j] = ccDict[g, cc]
end
##


scChar_ = []
for c in concepts
    sChar = Array{String}(undef, size(taxa_df, 1), length(sounds))
    for (i, l) in enumerate(taxa_df.Language), (j, s) in enumerate(sounds)
        g = taxa_df.Glottocode[i]
        sChar[i, j] = string(scDict[g, c, s])
    end
    push!(scChar_, sChar)
end

scMtx = hcat(scChar_...)

charMtx = hcat(ccMtx, scMtx)


##
taxa = string.(taxa_df.Language)


pad = maximum(length.(taxa)) + 5

nex = """
#Nexus

BEGIN DATA;
DIMENSIONS ntax=$(length(taxa)) NCHAR=$(size(charMtx, 2));
FORMAT DATATYPE=restriction GAP=? MISSING=- interleave=yes;
MATRIX

"""

for (i, l) in enumerate(taxa)
    global nex
    ln = rpad(l, pad) * join(charMtx[i, :]) * "\n"
    nex *= ln
end

nex *= """

;

END;
""";

##



open("../data/charMtx.nex", "w") do io
    write(io, nex)
end

##


constraints = []

for nd in glot.traverse()
    if !(nd.is_leaf() || nd.is_root())
        push!(
            constraints,
            join(nd.get_leaf_names(), " "),
        )
    end
end

##


constraints = []

for nd in glot.traverse()
    if !(nd.is_leaf() || nd.is_root())
        push!(
            constraints,
            join(nd.get_leaf_names(), " "),
        )
    end
end

##

mb = """
#Nexus
Begin MrBayes;
    execute ../data/charMtx.nex;
    charset cc = 1-$(size(ccMtx, 2));
    charset sc = $(size(ccMtx, 2)+1)-$(size(charMtx,2));
    partition dtype = 2:cc, sc;
    set partition = dtype;
    unlink Statefreq=(all) shape=(all);
    lset applyto=(all) rates=gamma;
    lset applyto=(1) coding=Noabsencesites;
    lset applyto=(2) coding=all;
"""

try
    mkdir("../output")
catch e
end

for (i, c) in enumerate(constraints)
    global mb
    ln = "    constraint c$i = $c;\n"
    mb *= ln
end

mb *= """
    prset topologypr = constraints($(join(["c$i" for i in 1:length(constraints)], ",")));
    prset brlenspr = clock:uniform;
    prset clockvarpr = igr;
    set beagleprecision=double;
    mcmcp stoprule=yes burninfrac=0.5 stopval=0.01 filename=../output/sample30 samplefreq=1000;
    mcmc ngen=100000000 nchains=2 nruns=2;
    sumt;
    q;
end;
""";

##


open("mbCommands.nex", "w") do io
    write(io, mb)
end
