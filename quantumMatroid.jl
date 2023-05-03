using Oscar: permutations

using Oscar
using Polymake
using Combinatorics

uni = uniform_matroid(2,4);
circ = [[1,4,5],[2,3,5],[1,2,3,4]];
ex = matroid_from_circuits(circ,collect(1:5));

matroid_groundset(ex)

function matroidIncidence(M::Matroid,t::String="circuits")
    rows = [];
    if (t === "circuits")
       rows = circuits(M);
    elseif (t === "hyperplanes")
        rows = hyperplanes(M);
    elseif (t === "flats")
        rows = flats(M);
    elseif (t === "bases")
        rows = bases(M);
    elseif (t === "independent_sets")
        rows = independent_sets(M);
    end

    groundset = matroid_groundset(M);
    incidence = zeros(Int64,length(rows),length(groundset));

    for i in 1:length(rows)
        for j in 1:length(groundset)
            if groundset[j] in rows[i]
                incidence[i,j] = 1;
            end
        end
    end
    
    return incidence;
end



function toAdjacency(Inc::Matrix{Int64})
    t = hcat(zeros(Int64,size(Inc,1),size(Inc,1)),Inc)
    t1 = hcat( transpose(Inc),zeros(Int64,size(Inc,2),size(Inc,2)))
    adj = vcat(t,t1)
    return adj
end

Inc = matroidIncidence(ex)

toAdjacency(Inc)






function autom_group(M::Matroid,t::String="circuits")
    if (t === "circuits")
        return automorphism_group(IncidenceMatrix(circuits(M)); action=:on_cols);
    elseif (t === "hyperplanes")
        return automorphism_group(IncidenceMatrix(hyperplanes(M)); action=:on_cols);
    elseif (t === "flats")
        fl = Vector{Vector{Int}}(flats(M))
        filter!(x->length(x) > 1,fl)
        return automorphism_group(IncidenceMatrix(fl); action=:on_cols);
    elseif (t === "bases")
        return automorphism_group(IncidenceMatrix(bases(M)); action=:on_cols);
    elseif (t === "independent_sets")
        ind = Vector{Vector{Int}}(independent_sets(M))
        filter!(x->length(x) > 1,ind)
        return automorphism_group(IncidenceMatrix(ind); action=:on_cols);
    else 
        println("Error: Invalid type");
    end
end

autom_group(ex,"circuits")
autom_group(uni,"bases")

graph =Graph{Undirected}(4);
edg = [[2,3],[1,2],[1,4],[3,4],[1,3]];
foreach(x->add_edge!(graph,x[1],x[2]),edg);

k4 = cycle_matroid(complete_graph(4));
bases(k4)
nonbases(k4)

b1 = bases(k4)[1]
nb1 = nonbases(k4)[1]

function matroidRelations(M::Matroid)
    nb = []
    b = []
    for base in bases(M)
        b = vcat(b,collect(permutations(base)))
    end
    b

    for nbase in nonbases(M)
        nb = vcat(nb,collect(permutations(nbase)))
    end
    nb
    r=rank(M)
    rels=[]
    for base in b
        for nbase in nb
            rel = []
            for i in 1:r
                push!(rel,(base[i],nbase[i]))
            end
            push!(rels,rel)
        end
    end
    for base in b
        for nbase in nb
            rel = []
            for i in 1:r
                push!(rel,(nbase[i],base[i]))
            end
            push!(rels,rel)
        end
    end
    return rels
end

matroidRelations(fano_matroid())
