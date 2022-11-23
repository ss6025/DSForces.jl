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


function ComputeDSForce(seq::Union{LongDNA{4}, LongSubSeq{DNAAlphabet{4}}})
    CSL = FindLongestComplementarySegment(seq)
    α = ComputeAlpha(seq)
    return CompSegmLenToDSForce(CSL, α, length(seq))
end

function ComputeDSForce(seq::AbstractString)
    # here I should check that no non-ACTG symbol is present
    DNA_seq = LongDNA{4}(seq)
    return ComputeDSForce(DNA_seq)
end

end
