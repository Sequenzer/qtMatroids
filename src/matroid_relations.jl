#!/usr/bin/julia
import AbstractAlgebra: QQ, elem_type
import AbstractAlgebra.Generic: FreeAssociativeAlgebra, FreeAssAlgElem, AhoCorasickAutomaton, insert_keyword!, normal_form

function add_new_relation!!(relations::Vector{AbstractAlgebra.Generic.FreeAssAlgElem{T}}, aut::AhoCorasickAutomaton, new_relation::FreeAssAlgElem{T}, relation_count=0::Int) where T
    normalized =new_relation# normal_form(new_relation, relations, aut)
    if !iszero(normalized)
        push!(relations, normalized)
        insert_keyword!(aut, normalized.exps[1], length(relations))
    end
    if relation_count % 100 == 0
#        interreduce!(relations)
    end
end

function matroid_relations(relation_indices::Vector{Any}, n::Int)
    generator_strings = String[]
    relation_count = 0
    for i in 1:n, j in 1:n
            push!(generator_strings, "u[$i,$j]")
    end
    A, g = FreeAssociativeAlgebra(QQ, generator_strings)
    u = Matrix{elem_type(A)}(undef, n, n)
    for i in 1:n, j in 1:n
            u[i, j] = g[(i-1)*n + j]
    end
    relations = elem_type(A)[]
    aut = AhoCorasickAutomaton(Vector{Int}[])
    for i in 1:n, j in 1:n
        new_relation = u[i, j] * u[i, j] - u[i, j]
        relation_count += 1
        if length(relations) == 0
            push!(relations, new_relation)
            insert_keyword!(aut, new_relation.exps[1], length(relations))
        else
            add_new_relation!!(relations, aut, new_relation, relation_count)
        end
            for k in 1:n
                    if k != j
                        relation_count += 1
                        new_relation = u[i,j] * u[i, k]
                        add_new_relation!!(relations, aut, new_relation, relation_count)
                        new_relation = u[j, i]*u[k, i]
                        relation_count += 1
                        add_new_relation!!(relations, aut, new_relation, relation_count)
#                                for l in 1:n
#                                        if adj[i, k] != adj[j, l]
#                                            new_relation = u[i, j]*u[k, l]
#                                            relation_count += 1
#                                            add_new_relation!!(relations, aut, new_relation, relation_count)
#
#                                        end
#                                end
                    end
            end
    end
    i = 1
    m = length(relation_indices)
    for relation in relation_indices
        if i % 5000 == 0
            percent = round(100*i/m,digits=2)
            println("i = $i, percent = $percent %")
        end
        i+=1
        temp = one(A)
        for gen in relation
            temp = temp * u[gen[1], gen[2]]
        end
        relation_count += 1
        add_new_relation!!(relations, aut, temp, relation_count)
    end
    return A, u, relations
    
end
function check_commutativity(u::Matrix{AbstractAlgebra.Generic.FreeAssAlgElem{Rational{BigInt}}}, gb::Vector{AbstractAlgebra.Generic.FreeAssAlgElem{Rational{BigInt}}}, aut::AhoCorasickAutomaton)
        for i in 1:size(u)[1]
                for j in 1:size(u)[2]
                        for k in 1:size(u)[1]
                                for l in 1:size(u)[2]
                                        if !iszero(normal_form(u[i, j]*u[k, l] - u[k, l] * u[i, j], gb, aut))
                                                return false
                                        end
                                end
                        end
                end
        end
        return true

end

function check_commutativity(u::Matrix{AbstractAlgebra.Generic.FreeAssAlgElem{Rational{BigInt}}}, gb::Vector{AbstractAlgebra.Generic.FreeAssAlgElem{Rational{BigInt}}})
        for i in 1:size(u)[1]
                for j in 1:size(u)[2]
                        for k in 1:size(u)[1]
                                for l in 1:size(u)[2]
                                        if !iszero(normal_form(u[i, j]*u[k, l] - u[k, l] * u[i, j], gb))
                                                return false
                                        end
                                end
                        end
                end
        end
        return true
end
