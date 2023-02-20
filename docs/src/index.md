```@meta
CurrentModule = DSForces
```

# DSForces

Documentation for [DSForces](https://github.com/adigioacchino/DSForces.jl), a Julia package to computed double stranded (DS) forces in nucleic-acid sequences.

## Installation
This package is not registered. Install with:

```julia
import Pkg
Pkg.add(url="https://github.com/adigioacchino/DSForces.jl")
```

## DS force definition
TBW (with the citation to the paper/bioRxiv).

## Examples
### Force computation for a given sequence
The only exported function is `ComputeDSForce`, that can be used to compute the DS force of any sequence, 
for instance the block
```julia
using DSForces
seq = "ACAACGTAACGGTCGAGTCG"
ComputeDSForce(seq)
```
will return a single `Float64` value that is the DS force associated to the sequence.
If the ranges of the two complementary (including Wobble pairs) sequences are needed, the parameter
`return_LCS_positions` must be set to `true`: 
```julia
using DSForces
seq = "ACAACGTAACGGTCGAGTCG"
ComputeDSForce(seq, return_LCS_positions=true)
```
The option `return_LCS_positions` will be used in all examples above, but can always be dropped if the
positions of the complementary sequences are not of interest.

### Force computation with constrained positions
Sometimes it can be useful to compute the DS force under the condition that one of the two subsequences 
forming the DS segment is in a given range.
This can be done with the argument `` such as in the following example
```julia
using DSForces
seq = "ACAACGTAACGGTCGAGTCG"
ComputeDSForce(seq, return_LCS_positions=true, seqA_range=1:10)
```
Here we are restricting the position of one of the two subsequences to the range `1:10`, and this 
changes the length of the longest complementary segment (now 2, compare with 3 in the unrestricted case), 
as well as the value of the DS force. 
Notice that a different equation is used to compute the DS force from the length of the longest DS segment
if `seqA_range` is specified, to account for the new constraint.

### Force computation with sliding windows
For very long sequences it could be interesting to compute the DS force for a window of fixed lenght that
slides along the sequence. 
The following code block demonstrates the usage of the argument `sliding_window_length` that allows for this:
```julia
using DSForces
pre_seq = "ACAACGTAACGGTCGAGTCG"
seq = pre_seq^2
ComputeDSForce(seq, return_LCS_positions=true, sliding_window_length=length(pre_seq))
```
Here the length of the sliding window is the length of `pre_seq`, which is identical to the examples above.
Therefore the first element of the 3 vectors that the function outputs will be identical to the output of the 
second example in section [Force computation for a given sequence](@ref), and all other values of these vectors 
will be the same quantities computed for the other sliding windows of the same length along the sequence. 
Notice that the stride is fixed to 1.

## API

```@index
```

```@autodocs
Modules = [DSForces]
```
