VERSION < v"0.7.0-beta2.199" && __precompile__()

module MDCT
using Compat
export mdct, imdct

if VERSION < v"0.7.0-DEV.602"
    using Base.FFTW
    import Base.FFTW: fftwNumber, r2r, r2r!, REDFT11
else
    using FFTW
    import FFTW: fftwNumber, r2r, r2r!, REDFT11
end

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

end # MDCT
