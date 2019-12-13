const N = 624
const M = 397
const UPPER_MASK = 0x80000000
const LOWER_MASK = 0x7fffffff

mutable struct MT19937
    mt::Vector{UInt32}
    mti::Int
    function MT19937(x::Vector{UInt32}, i::Int)
        @assert length(x) == N
        new(x, i)
    end
end

@inline mt_magic(y) = ((y % Int32) << 31 >> 31) & 0x9908b0df
"Get a random `UInt32` number from a `MT19937` object."
@inline function mt_get(r::MT19937)
    mt = r.mt
    if r.mti > N
        @inbounds for i in 1:N-M
            y = (mt[i] & UPPER_MASK) | (mt[i+1] & LOWER_MASK)
            mt[i] = mt[i + M] ⊻ (y >> 1) ⊻ mt_magic(y)
        end
        @inbounds for i in N-M+1:N-1
            y = (mt[i] & UPPER_MASK) | (mt[i+1] & LOWER_MASK)
            mt[i] = mt[i + M - N] ⊻ (y >> 1) ⊻ mt_magic(y)
        end
        @inbounds begin
            y = (mt[N] & UPPER_MASK) | (mt[1] & LOWER_MASK)
            mt[N] = mt[M] ⊻ (y >> 1) ⊻ mt_magic(y)
        end
        r.mti = 1
    end
    k = mt[r.mti]
    k ⊻= (k >> 11)
    k ⊻= (k << 7) & 0x9d2c5680
    k ⊻= (k << 15) & 0xefc60000
    k ⊻= (k >> 18)

    r.mti += 1
    k
end
