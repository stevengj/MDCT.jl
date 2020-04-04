VERSION < v"0.7.0-beta2.199" && __precompile__()

module MDCT

using LinearAlgebra
import AbstractFFTs
export mdct, imdct, plan_mdct, plan_imdct
import Base: *, size
import LinearAlgebra: mul!

using FFTW
import FFTW: fftwNumber, r2r, r2r!, REDFT11, plan_r2r!, plan_r2r

fftwsimilar(X::AbstractArray{T}, sz) where {T<:fftwNumber} = Array{T}(undef, sz...)
fftwsimilar(X::AbstractArray{T}, sz) where {T<:Real} = Array{Float64}(undef, sz...)
fftwsimilar(X::AbstractArray{T}, sz) where {T<:Complex} = Array{Complex128}(undef, sz...)

# The following two routines compute the MDCT and IMDCT via
# a type-IV DCT (FFTW's REDFT11 r2r transform).  For a review
# of this relationship, see the notes I posted on Wikipedia:
#     http://en.wikipedia.org/wiki/Modified_discrete_cosine_transform

function mdct(X::AbstractVector{T}) where {T<:Number}
    sz = length(X)
    if isodd(sz)
        throw(ArgumentError("mdct requires an even-length vector"))
    end
    N = div(sz, 2);
    Y = fftwsimilar(X, N)
    if isodd(N)
        throw(ArgumentError("mdct requires a multiple-of-4 vector length"))
        # FIXME: handle odd case via DCT-III?
    else
        N2 = div(N,2)
        for i = 1:N2
            Y[i] = -0.5 * (X[(3*N2+1)-i] + X[3*N2+i])
            Y[N2+i] = -0.5 * (X[(2*N2+1)-i] - X[i])
        end
        r2r!(Y, REDFT11)
    end
    return Y
end

function imdct(X::StridedVector{T}) where {T<:Number}
    N = length(X)
    Z = fftwsimilar(X, 2*N)
    if isodd(N)
        throw(ArgumentError("imdct requires an even vector length"))
        # FIXME: handle odd case via DCT-II?
    else
        Y = r2r(X, REDFT11)
        N2 = div(N,2)
        s = 0.5 / N
        for i = 1:N2
            Z[N+1-i] = -(Z[i] = Y[N2+i]*s)
            Z[N+N2+1-i] = (Z[N+N2+i] = -Y[i]*s)
        end
    end
    return Z
end

imdct(X::AbstractVector{T}) where {T<:Number} =
    imdct(copy!(fftwsimilar(X, size(X)), X))

mutable struct Plan{T<:fftwNumber, N, inv} <: AbstractFFTs.Plan{T}
    plan::FFTW.r2rFFTWPlan{T} # plan_r2r! REDFT11 plan
    pinv::AbstractFFTs.Plan{T}
    Plan{T, N, inv}(plan) where {T, N, inv} = new(plan)
end

size(p::Plan{T, N, true}) where {T<:Number, N} = (div(N,2), N)
size(p::Plan{T, N, false}) where {T<:Number, N} = (2*N, N)

function mul!(Y::StridedVector{T}, p::Plan{T, N, false}, X::AbstractVector{T}) where {T<:Number, N}
    @boundscheck (length(X) == 2*N && length(Y) == N) || throw(DimensionMismatch())
    N2 = div(N,2)
    @inbounds for i = 1:N2
        Y[i] = -0.5 * (X[(3*N2+1)-i] + X[3*N2+i])
        Y[N2+i] = -0.5 * (X[(2*N2+1)-i] - X[i])
    end
    return p.plan * Y
end

function mul!(Z::StridedVector{T}, p::Plan{T, N, true}, X::StridedVector{T}) where {T<:Number, N}
    @boundscheck length(Z) == N || throw(DimensionMismatch())
    Y = p.plan*X
    N1 = div(N,2)
    N2 = div(N1,2)
    s = 0.5 / N1
    @inbounds for i = 1:N2
        Z[N1+1-i] = -(Z[i] = Y[N2+i]*s)
        Z[N1+N2+1-i] = (Z[N1+N2+i] = -Y[i]*s)
    end
    return Z
end

function *(p::Plan{T, N}, X::AbstractVector{T}) where {T<:Number, N}
    Y = fftwsimilar(X, N)
    return mul!(Y, p, X)
end

function plan_mdct(X::AbstractVector{T}, flags::Integer=FFTW.ESTIMATE, timelimit::Real=FFTW.NO_TIMELIMIT) where {T<:Number}
    sz = length(X)
    if mod(sz, 4) != 0
        throw(ArgumentError("mdct requires a multiple-of-4 vector length"))
        # FIXME: handle odd case via DCT-III?
    end
    N = div(sz, 2);
    return Plan{T, N, false}(plan_r2r!(fftwsimilar(X, N), REDFT11; flags=flags, timelimit=timelimit))
end

function plan_imdct(X::AbstractVector{T}, flags::Integer=FFTW.ESTIMATE, timelimit::Real=FFTW.NO_TIMELIMIT) where {T<:Number}
    sz = length(X)
    if isodd(sz)
        throw(ArgumentError("imdct requires an even vector length"))
        # FIXME: handle odd case via DCT-II?
    end
    return Plan{T, 2sz, true}(plan_r2r(fftwsimilar(X, sz), REDFT11; flags=flags, timelimit=timelimit))
end

end # MDCT
