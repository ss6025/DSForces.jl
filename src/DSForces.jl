module DSForces

using BioSequences
using StatsBase
using Random

export ComputeDSForce

const allowed_alphabet = [DNA_A, DNA_C, DNA_G, DNA_T]
const allowed_pairs = [[DNA_A, DNA_T], [DNA_C, DNA_G], [DNA_G, DNA_T], [DNA_T, DNA_A], [DNA_G, DNA_C], [DNA_T, DNA_G]]


"""
    CompSegmLenToDSForce(CSL::Int, alpha::AbstractFloat, seqlen::Int, 
                              seqlen_forced::Union{Int, Missing}=seqlen, c0::AbstractFloat=-2.2)

Given the length of the longest complementary segmenet `CSL`, the value of the probability
`alpha` that two randomly chosen nucleotide in the sequence under analysis form an 
allowed pair, and the length of the sequence under analysis `lenseq`, it computes the DS force.
The value of the parameter `c0` has been fixed by synthetic tests, see the accompanying paper.
If one of the two complementary sequences is forced to be within an interval (whose length is 
`seqlen_forced`), the DS force computation is adapted to this case. 
"""
function CompSegmLenToDSForce(CSL::Int, alpha::AbstractFloat, seqlen::Int, 
                              seqlen_forced::Union{Int, Missing}=seqlen, c0::AbstractFloat=-2.2)
    if CSL < 2
        return NaN
    else
        return log(1 / alpha) - (log(seqlen) + log(seqlen_forced)) / (CSL - c0)
    end
end


"""
    ComputeAlpha(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}})

Compute the probability ("α") that two randomly chosen nucleotide in `seq` form an 
allowed pair (defined by the global variable `allowed_pairs`).
"""
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


"""
    FindLongestComplementarySegment(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}}, 
                                         return_coords::Bool=false; rand_choice::Bool=false, 
                                         seed::Union{Missing, Int}=missing)
                                         return_coords::Bool=false)

Find the two subsequences of `seq` that can 
form the longest possible double-strand segement (including Watson-Crick and Wobble pairs).
If `return_coords` the ranges encompassing the two 'complementary' sequences are returned,
otherwise only the length of the segment is returned.
When different ranges are possible, those closest to the end of the sequence ("downstream") 
are chosen, unless `rand_choice` is set to `true`, in which case a random choice is made. 
"""
function FindLongestComplementarySegment(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}}, 
                                         return_coords::Bool=false; rand_choice::Bool=false, 
                                         seed::Union{Missing, Int}=missing)
    pattern_seq = LongDNA{4}(replace(String(seq), "A" => "T", "C" => "G", "G" => "Y", "T"=>"R"))
    L = length(seq)
    maxL = 2
    start_segm_poss = [1]
    @inbounds for i in 1:L
        for k in i+maxL-1:L
            if !(occursin(ExactSearchQuery(reverse(pattern_seq[i:k]), iscompatible), @view seq[k+1:end]))
                if (k - i) > maxL
                    maxL = k - i
                    start_segm_poss = [i]
                elseif (k - i) == maxL 
                    push!(start_segm_poss, i)
                end
                break
            end
        end
    end
    if return_coords
        if rand_choice
            ismissing(seed) ? rngX = Xoshiro() : rngX = Xoshiro(seed)
            start_segm_pos = rand(rngX, start_segm_poss) # random choice of starting position
        else
            start_segm_pos = start_segm_poss[end] # downstream choice of starting position
        end
        eqs = ExactSearchQuery(reverse(pattern_seq[start_segm_pos:start_segm_pos+maxL-1]), iscompatible)
        downstream_poss = findall(eqs, seq[start_segm_pos+maxL:end]) 
        if length(downstream_poss) == 0
            return L:L, L:L
        else
            if rand_choice
                return start_segm_pos:start_segm_pos+maxL-1, rand(rngX, downstream_poss) .+ (start_segm_pos+maxL-1) # random choice of downstream position
            else
                return start_segm_pos:start_segm_pos+maxL-1, downstream_poss[1] .+ (start_segm_pos+maxL-1) # closest choice of downstream position
            end
        end
    else
        return maxL
    end
end


"""
    FindLongestComplementarySegment(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}}, 
                                         seqA_range::UnitRange{Int}, return_coords::Bool=false)

Find the longest subsequence of `seq` that is within `seqA_range` and can form a double-strand 
segement (including Watson-Crick and Wobble pairs) with another subsequence of `seq`.
If `return_coords` the ranges encompassing the two 'complementary' sequences are returned,
otherwise only the length of the segment is returned.
"""
function FindLongestComplementarySegment(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}}, 
                                         seqA_range::UnitRange{Int}, return_coords::Bool=false)
    pattern_seq = LongDNA{4}(replace(String(seq), "A" => "T", "C" => "G", "G" => "Y", "T"=>"R"))
    L = length(seq)
    maxL = 2
    start_segm_pos = seqA_range[1]
    @inbounds for i in seqA_range
        for k in i+maxL-1:seqA_range[end]
            esq = ExactSearchQuery(reverse(pattern_seq[i:k]), iscompatible)
            c1 = occursin(esq, @view seq[1:i-1]) # paired segment is before seqA
            c2 = occursin(esq, @view seq[k+1:end]) # paired segment is after seqA
            if (!(c1) & !(c2))
                if (k - i) > maxL
                    maxL = k - i
                    start_segm_pos = i
                elseif (k - i) == maxL 
                    start_segm_pos = i
                end
                break
            elseif k == seqA_range[end]
                if (k - i + 1) > maxL
                    maxL = k - i + 1
                    start_segm_pos = i
                elseif (k - i + 1) == maxL 
                    start_segm_pos = i
                end
            end
        end
    end
    if return_coords
        eqs = ExactSearchQuery(reverse(pattern_seq[start_segm_pos:start_segm_pos+maxL-1]), iscompatible)
        pre_match = findlast(eqs, seq[1:start_segm_pos-maxL])
        post_match = findfirst(eqs, seq[start_segm_pos+maxL:end])
        if isnothing(pre_match) & isnothing(post_match)
            return L:L, L:L
        elseif !(isnothing(post_match))
            return start_segm_pos:start_segm_pos+maxL-1, post_match .+ (start_segm_pos+maxL-1)
        else
            return pre_match, start_segm_pos:start_segm_pos+maxL-1
        end
    else
        return maxL
    end
end


"""
    FindLongestComplementarySegmentLast(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}})

Find the longest subsequence of `seq` that ends at the end of `seq` and can form a double-strand 
segement (including Watson-Crick and Wobble pairs) with another subsequence of `seq`, and return
the coordinates of these two subsequences.
NOTE: This function is only used to speed up the sliding computation of DS forces.
"""
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


"""
    FindLongestComplementarySegmentSliding(seq::LongDNA{4}, 
    sliding_window_length::Int=min(length(seq), 3000))

Does the same as `FindLongestComplementarySegment` on sliding windows of length
`sliding_window_length` and with stride 1. It returns four vectors, each with one 
value for each sliding window: 
1, the lengths of the longest complementary segment; 
2, the values of α (see `ComputeAlpha`); 
3, the position of the upstream complementary subsequence;
4, the position of the downstream complementary subsequence.
"""
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


"""
    ComputeDSForce(
            seq::Union{AbstractString, LongDNA{4}, LongSubSeq{DNAAlphabet{4}}};
            return_LCS_positions::Bool=false, 
            sliding_window_length::Union{Int, Missing}=missing,
            seqA_range::Union{UnitRange{Int}, Missing}=missing)

Given the sequence `seq`, compute its DS force.
If `return_LCS_positions`, the positions of the subsequences forming the 
longest double-strand segment are returned.
If `sliding_window_length` is provided, the DS force computation is done for 
sliding windows of that length, with stride 1.
If a range is specified as `seqA_range`, one of the two fully complementary 
is forced to be within the provided range.

WARNGING: specifying a `seqA_range` has no effect if also a 
`sliding_window_length` is provided.
"""
function ComputeDSForce(
            seq::Union{AbstractString, LongDNA{4}, LongSubSeq{DNAAlphabet{4}}};
            return_LCS_positions::Bool=false, 
            sliding_window_length::Union{Int, Missing}=missing,
            seqA_range::Union{UnitRange{Int}, Missing}=missing)
    if typeof(seq) <: AbstractString
        return ComputeDSForce(LongDNA{4}(seq); return_LCS_positions=return_LCS_positions, 
                sliding_window_length=sliding_window_length, 
                seqA_range=seqA_range)
    else
        # return missing if ambiguous nucletides are present
        (sum(isambiguous.(seq)) > 0) && return missing
        if ismissing(sliding_window_length)
            α = ComputeAlpha(seq)
            if return_LCS_positions
                if ismissing(seqA_range)
                    rangeA, rangeB = FindLongestComplementarySegment(seq, true)
                    CSL = length(rangeA)
                    return CompSegmLenToDSForce(CSL, α, length(seq)), rangeA, rangeB
                else
                    rangeA, rangeB = FindLongestComplementarySegment(seq, seqA_range, true)
                    CSL = length(rangeA)
                    return CompSegmLenToDSForce(CSL, α, length(seq), length(seqA_range)), rangeA, rangeB
                end
            else
                if ismissing(seqA_range)
                    CSL = FindLongestComplementarySegment(seq)
                    return CompSegmLenToDSForce(CSL, α, length(seq))
                else
                    CSL = FindLongestComplementarySegment(seq, seqA_range)
                    return CompSegmLenToDSForce(CSL, α, length(seq), length(seqA_range))
                end
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
end


end
