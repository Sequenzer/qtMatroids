using Oscar
using Polymake
using Combinatorics

include("./matroid_relations.jl")
include("./utils.jl")



struct MultiSetMatroid
    classic::Matroid
end

Oscar.bases(M::MultiSetMatroid) = bases(M.classic)
Oscar.independent_sets(M::MultiSetMatroid)= independent_sets(M.classic)

function Oscar.circuits(M::MultiSetMatroid)
    matroid = M.classic
    groundset = matroid_groundset(matroid)
    cc = Oscar.circuits(matroid)
    lops = Oscar.loops(matroid)
    for ele in 1:length(matroid)
        ele in lops && continue
        for ele2 in ele:length(matroid)
            ele2 in lops && continue
            [ele,ele2] in cc && continue
            push!(cc,[ele,ele2])
        end
    end
    return cc
end

A = [1,2]
unique(sort.([[1,2],[2,1]]))

unique(collect(powerset([1,2,1,2])))
reduce(vcat,[A for _ in 1:length(A)])

function matroidRelations(M::MultiSetMatroid,structure::Symbol=:bases)

    matroid = M.classic
    n = length(matroid)
    grdSet = matroid_groundset(matroid)

    b =  [[] for _ in 1:n]
    nb =  [[] for _ in 1:n]
    rels=[]

    sets  = eval(structure)(M)
    sizes = unique(map(x -> length(x), sets))

    for size in sizes 
        tempGrdSet = reduce(vcat,[grdSet for i in 1:n])
        powerSet = unique(sort.(powerset(tempGrdSet,size,size)))

        setsOfSize = filter(x->length(x)==size,sets)
        nonSets = setdiff(powerSet,setsOfSize) 

        for set in setsOfSize
            b[size] = vcat(b[size],collect(permutations(set)))
        end
        for nonset in nonSets
            nb[size] = vcat(nb[size],collect(permutations(nonset)))
        end
        for set in b[size]
            for nonset in nb[size]
                rel = []
                for i in 1:size
                    push!(rel,(set[i],nonset[i]))
                end
                push!(rels,rel)
                rel = []
                for i in 1:size
                    push!(rel,(nonset[i],set[i]))
                end
                push!(rels,rel)
            end
        end

    end

    return unique(rels)
end

#=
B = [[1],[2]]
n_Elements = 5 
ex = matroid_from_bases(B,n_Elements)    
multiEx = MultiSetMatroid(ex) 
circuits(multiEx)

relations_b = matroidRelations(multiEx,:bases)
result_b = matroid_relations(relations_b,n_Elements)
result[3]


relations_c = matroidRelations(multiEx,:circuits)
result_c = matroid_relations(relations_c,n_Elements)
result_c[3]

=#


#=

multiFano = MultiSetMatroid(fano_matroid())
n_Elements = length(fano_matroid())
matroidRelations(multiFano,:bases)
circuits(multiFano)
matroidRelations(multiFano,:circuits)


relations_b = matroidRelations(multiFano,:bases)
result_b = matroid_relations(relations_b,n_Elements)
result[3]


relations_c = matroidRelations(multiFano,:circuits)
result_c = matroid_relations(relations_c,n_Elements)
result_c[3]


=#









#=
B = [[1],[2]]
n_Elements = 5 
ex = matroid_from_bases(B,n_Elements)
circuits(ex)

#Get the Adjacencygraph for the bipartite Graph
Inc = matroidIncidence(ex)
toAdjacency(Inc)
#Get the Relations definign Aut^\+(M)

=#


