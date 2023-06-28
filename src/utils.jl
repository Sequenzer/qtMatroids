using Oscar
using Polymake
using Combinatorics


function matroidIncidence(M::Matroid,t::Symbol=:circuits)
    rows = eval(t)(M);

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


