module MDCT
export mdct, imdct

import Base.FFTW.fftwNumber, Base.FFTW.r2r, Base.FFTW.r2r!, Base.FFTW.REDFT11

fftwsimilar{T<:fftwNumber}(X::AbstractArray{T}, sz) = Array(T, sz...)
fftwsimilar{T<:Real}(X::AbstractArray{T}, sz) = Array(Float64, sz...)
fftwsimilar{T<:Complex}(X::AbstractArray{T}, sz) = Array(Complex128, sz...)

# The following two routines compute the MDCT and IMDCT via
# a type-IV DCT (FFTW's REDFT11 r2r transform).  For a review
# of this relationship, see the notes I posted on Wikipedia:
#     http://en.wikipedia.org/wiki/Modified_discrete_cosine_transform

function mdct{T<:Number}(X::AbstractVector{T})
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

function imdct{T<:Number}(X::StridedVector{T})
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

imdct{T<:Number}(X::AbstractVector{T}) = imdct(copy!(fftwsimilar(X, size(X)),
                                                     X))

end # MDCT
