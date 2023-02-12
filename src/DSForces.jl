module DSForces

using BioSequences
using StatsBase

export ComputeDSForce

const allowed_alphabet = [DNA_A, DNA_C, DNA_G, DNA_T]
const allowed_pairs = [[DNA_A, DNA_T], [DNA_C, DNA_G], [DNA_G, DNA_T], [DNA_T, DNA_A], [DNA_G, DNA_C], [DNA_T, DNA_G]]


function CompSegmLenToDSForce(CSL::Int, alpha::AbstractFloat, seqlen::Int, c0::AbstractFloat=-2.2)
    if CSL < 2
        return NaN
    else
        return log(1 / alpha) - 2 * log(seqlen) / (CSL - c0)
    end
end


function ComputeAlpha(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}})
    seqlen = length(seq)
    nt_usage = countmap(seq)
    ks_ntu = keys(nt_usage)
    if sum([nt in ks_ntu ? nt_usage[nt] : 0 for nt in allowed_alphabet]) != seqlen
        return 0.
    else
        alpha = sum([((p[1] in ks_ntu) & (p[2] in ks_ntu)) ? nt_usage[p[1]] * nt_usage[p[2]] / seqlen^2 : 0 for p in allowed_pairs])
        return alpha
    end
end


function FindLongestComplementarySegment(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}}, 
                                         return_coords::Bool=false)
    pattern_seq = LongDNA{4}(replace(String(seq), "A" => "T", "C" => "G", "G" => "Y", "T"=>"R"))
    L = length(seq)
    maxL = 2
    start_segm_pos = 1
    @inbounds for i in 1:L
        for k in i+maxL-1:L
            if !(occursin(ExactSearchQuery(reverse(pattern_seq[i:k]), iscompatible), @view seq[k+1:end]))
                if (k - i) > maxL
                    maxL = k - i
                    start_segm_pos = i
                elseif (k - i) == maxL 
                    start_segm_pos = i
                end
                break
            end
        end
    end
    if return_coords
        eqs = ExactSearchQuery(reverse(pattern_seq[start_segm_pos:start_segm_pos+maxL-1]), iscompatible)
        if isnothing(findfirst(eqs, seq[start_segm_pos+maxL:end]))
            return L:L, L:L
        else
            return start_segm_pos:start_segm_pos+maxL-1, findfirst(eqs, seq[start_segm_pos+maxL:end]) .+ (start_segm_pos+maxL-1)
        end
    else
        return maxL
    end
end

function FindLongestComplementarySegmentLast(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}})
    pattern_seq = LongDNA{4}(replace(String(seq), "A" => "T", "C" => "G", "G" => "Y", "T"=>"R"))
    L = length(seq)
    maxL = 2
    start_segm_pos = L - 2
    @inbounds for k in maxL:L
        if !(occursin(ExactSearchQuery(reverse(pattern_seq[L-k:L]), iscompatible), @view seq[1:L-k-1]))
            if k+1 >= maxL
                maxL = k + 1
                start_segm_pos = L-k
            end
            break
        end
    end
    eqs = ExactSearchQuery(reverse(pattern_seq[start_segm_pos+1:L]), iscompatible)
    if isnothing(findfirst(eqs, seq[1:start_segm_pos]))
        return L:L, L:L
    else
        return findfirst(eqs, seq[1:start_segm_pos]), start_segm_pos+1:L
    end
end

function FindLongestComplementarySegmentSliding(seq::LongDNA{4}, 
    sliding_window_length::Int=min(length(seq), 3000))
    LCSLs = Int[]
    alphas = Float64[]
    range1s = UnitRange{Int64}[]
    range2s = UnitRange{Int64}[]
    global_pos = 1
    while global_pos <= length(seq)-sliding_window_length+1
        sliding_seq = seq[global_pos:global_pos+sliding_window_length-1]        
        nt_usage = countmap(sliding_seq)
        ks_ntu = keys(nt_usage)
        if sum([nt in ks_ntu ? nt_usage[nt] : 0 for nt in allowed_alphabet]) != sliding_window_length
            push!(LCSLs, 0)
            push!(alphas, 0)
            push!(range1s, 1:1)
            push!(range2s, 1:1)
            global_pos += 1
            continue
        end
        t_alpha = sum([((p[1] in ks_ntu) & (p[2] in ks_ntu)) ? nt_usage[p[1]] * nt_usage[p[2]] / sliding_window_length^2 : 0 for p in allowed_pairs])
        push!(alphas, t_alpha)
        if global_pos == 1
            t_range1, t_range2 = FindLongestComplementarySegment(sliding_seq, true)
            push!(range1s, t_range1 .+ (global_pos-1))
            push!(range2s, t_range2 .+ (global_pos-1))
        else
            pre_range1 = range1s[end]
            if pre_range1[1] < global_pos
                t_range1, t_range2 = FindLongestComplementarySegment(sliding_seq, true)
                push!(range1s, t_range1 .+ (global_pos-1))
                push!(range2s, t_range2 .+ (global_pos-1))
            else
                pre_length = LCSLs[end]
                lastpos_range1, lastpos_range2 = FindLongestComplementarySegmentLast(sliding_seq)
                if length(lastpos_range1) > pre_length
                    t_range1, t_range2 = lastpos_range1, lastpos_range2
                    push!(range1s, lastpos_range1 .+ (global_pos-1))
                    push!(range2s, lastpos_range2 .+ (global_pos-1))
                else
                    push!(range1s, range1s[end])
                    push!(range2s, range2s[end])
                end
            end
        end
        push!(LCSLs, length(range1s[end]))
        global_pos += 1
    end    
    return LCSLs, alphas, range1s, range2s
end


function ComputeDSForce(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}};
                        return_LCS_positions::Bool=false, 
                        sliding_window_length::Union{Int, Missing}=missing)
    if ismissing(sliding_window_length)
        α = ComputeAlpha(seq)
        if return_LCS_positions
            rangeA, rangeB = FindLongestComplementarySegment(seq, true)
            CSL = length(rangeA)
            return CompSegmLenToDSForce(CSL, α, length(seq)), rangeA, rangeB
        else
            CSL = FindLongestComplementarySegment(seq)
            return CompSegmLenToDSForce(CSL, α, length(seq))
        end
    else
        CSLs, alphas, range1s, range2s = FindLongestComplementarySegmentSliding(seq, sliding_window_length)
        DSFs = [CompSegmLenToDSForce(x, y, sliding_window_length) for (x,y) in zip(CSLs, alphas)]
        if return_LCS_positions
            return DSFs, range1s, range2s
        else
            return DSFs
        end
    end
end

function ComputeDSForce(seq::AbstractString;
                        return_LCS_positions::Bool=false, 
                        sliding_window_length::Union{Int, Missing}=missing)
    # here I should check that no non-ACTG symbol is present
    DNA_seq = LongDNA{4}(seq)
    return ComputeDSForce(DNA_seq; return_LCS_positions=return_LCS_positions, 
                            sliding_window_length=sliding_window_length)
end


end
