
"""
    coordinates(x)

Return the spatial coordinates of `x`.
"""
coordinates(x::AbstractArray) = Iterators.product(spatial_indices(x))

# TODO document
""" AbstractCoordinateSystem """
abstract type AbstractCoordinateSystem end

struct XYZCoordinateSystem <: AbstractCoordinateSystem end

@inline function orientation_name(s::AbstractCoordinateSystem, x::Tuple{Any,Any})
    (orientation_name(s, getfield(x, 1, false)),
     orientation_name(s, getfield(x, 2, false)))
end
@inline function orientation_name(s::AbstractCoordinateSystem, x::Tuple{Any,Any,Any})
    (orientation_name(s, getfield(x, 1, false)),
     orientation_name(s, getfield(x, 2, false)),
     orientation_name(s, getfield(x, 3, false)))
end

@inline function orientation_name(::XYZCoordinateSystem, x::Int)
    if x === 1
        return :left_to_right
    elseif x === -1
        return :right_to_left
    elseif x === 2
        return :back_to_front
    elseif x === -2
        return :front_to_back
    elseif x === 3
        return :bottom_to_top
    elseif x === -3
        return :top_to_bottom
    else
        error("$x does not map to a dimension name.")
    end
end

"""
    orientation(x)

Returns a tuple providing the orientation of `x`.
"""
@inline orientation(x) = _orientation(spatial_order(x), x)
_orientation(dims::Tuple, x) = dims
# if all dimnames are `:_` then there aren't any actual names and we need to derive from
# spatial_directions
@inline function _orientation(dims::Tuple{Vararg{StaticSymbol{:_}}}, x)
    orientation_name(AbstractCoordinateSystem(x), _dir2ori(spatial_directions(x)))
end
@inline function _dir2ori(x::Tuple{Tuple{Any,Any},Tuple{Any,Any}})
    x1 = getfield(x, 1, false)
    x2 = getfield(x, 2, false)
    return _dir2ori(
        getfield(x1, 1, false), getfield(x1, 2, false), 0,
        getfield(x2, 1, false), getfield(x2, 2, false), 0,
        0, 0, 0
    )
end
@inline function _dir2ori(x::Tuple{Tuple{Any,Any,Any},Tuple{Any,Any,Any},Tuple{Any,Any,Any}})
    x1 = getfield(x, 1, false)
    x2 = getfield(x, 2, false)
    x3 = getfield(x, 3, false)
    return _dir2ori(
        getfield(x1, 1, false), getfield(x1, 2, false), getfield(x1, 3, false),
        getfield(x2, 1, false), getfield(x2, 2, false), getfield(x2, 3, false),
        getfield(x3, 1, false), getfield(x3, 2, false), getfield(x3, 3, false),
    )
end


@inline function _det(r11, r12, r13, r21, r22, r23, r31, r32, r33)
    r11 * r22 * r33 -
    r11 * r32 * r23 -
    r21 * r12 * r33 +
    r21 * r32 * r13 +
    r31 * r12 * r23 -
    r31 * r22 * r13
end

@inline function _mul_trace(
    x11, x12, x13, x21, x22, x23, x31, x32, x33,
    y11, y12, y13, y21, y22, y23, y31, y32, y33
)

    return (x11 * y11 + x12 * y21 + x13 * y31) +  # z11
           (x21 * y12 + x22 * y22 + x23 * y32) +  # z22
           (x31 * y13 + x32 * y23 + x33 * y33)    # z33
end

function _dir2ori(xi, xj, xk, yi, yj, yk, zi, zj, zk)
    # Normalize column vectors to get unit vectors along each ijk-axis
    # normalize i axis
    val = sqrt(xi*xi + yi*yi + zi*zi)
    if val == 0
        error("Invalid rotation directions.")
    end
    xi /= val
    yi /= val
    zi /= val

    # normalize j axis
    val = sqrt(xj*xj + yj*yj + zj*zj)
    if val == 0
        error("Invalid rotation directions.")
    end
    xj /= val
    yj /= val
    zj /= val

    # orthogonalize j axis to i axis, if needed
    val = xi*xj + yi*yj + zi* zj  # dot product between i and j
    if abs(val) > .0001
        xj -= val*xi
        yj -= val*yi
        zj -= val*zi

        val = sqrt(xj*xj + yj*yj + zj*zj)  # must renormalize
        if val == 0
            error("The first and second dimensions cannot be parallel.")
        end
        xj /= val
        yj /= val
        zj /= val
    end

    # normalize k axis; if it is zero, make it the cross product i x j
    val = sqrt(xk*xk + yk*yk + zk*zk)
    if val == 0
        xk = yi*zj-zi*yj
        yk = zi*xj-zj*xi
        zk = xi*yj-yi*xj
    else
        xk = xk/val
        yk = yk/val
        zk = zk/val
    end

    # orthogonalize k to i
    val = xi*xk + yi*yk + zi*zk  # dot product between i and k
    if abs(val) > 0.0001
        xk -= val*xi
        yk -= val*yi
        zk -= val*zi

        # must renormalize
        val = sqrt(xk*xk + yk*yk + zk*zk)
        if val == 0
            return 0  # I think this is suppose to be an error output
        end
        xk /= val
        yk /= val
        zk /= val
    end

    # orthogonalize k to j */
    val = xj*xk + yj*yk + zj*zk  # dot product between j and k
    if abs(val) > 0.0001
        xk -= val*xj
        yk -= val*yj
        zk -= val*zj

        val = sqrt(xk*xk + yk*yk + zk*zk)
        if val == 0
            return 0  # bad
        end
        xk /= val
        yk /= val
        zk /= val
    end

    # at this point Q is the rotation matrix from the (i,j,k) to (x,y,z) axes
    detQ = _det(xi, xj, xk, yi, yj, yk, zi, zj, zk)
    # if( detQ == 0.0 ) return ; /* shouldn't happen unless user is a DUFIS */

    # Build and test all possible +1/-1 coordinate permutation matrices P;
    # then find the P such that the rotation matrix M=PQ is closest to the
    # identity, in the sense of M having the smallest total rotation angle.

    # Despite the formidable looking 6 nested loops, there are
    # only 3*3*3*2*2*2 = 216 passes, which will run very quickly.
    vbest = -666
    ibest = pbest=qbest=rbest= 1
    jbest = 2
    kbest = 3
    for (i, j, k) in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
        for p in (-1, 1)           # p,q,r are -1 or +1
            for q in (-1, 1)       # and go into rows 1,2,3
                for r in (-1, 1)
                    p11, p12, p13 = _nval_other_zero(i, p)
                    p21, p22, p23 = _nval_other_zero(j, q)
                    p31, p32, p33 = _nval_other_zero(k, r)
                    #=
                    P[1,i] = p
                    P[2,j] = q
                    P[3,k] = r
                    detP = det(P)  # sign of permutation
                    =#
                    detP = _det(p11, p12, p13, p21, p22, p23, p31, p32, p33)
                    # doesn't match sign of Q
                    if detP * detQ >= 0.0
                        # angle of M rotation = 2.0 * acos(0.5 * sqrt(1.0 + trace(M)))
                        # we want largest trace(M) == smallest angle == M nearest to I
                        val = _mul_trace(
                            p11, p12, p13, p21, p22, p23, p31, p32, p33,
                            xi, xj, xk, yi, yj, yk, zi, zj, zk
                        )
                        if val > vbest
                            vbest = val
                            ibest = i
                            jbest = j
                            kbest = k
                            pbest = p
                            qbest = q
                            rbest = r
                        end
                    end
                end
            end
        end
    end
    # At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.

    # The matrix P that corresponds is the best permutation approximation
    # to Q-inverse; that is, P (approximately) takes (x,y,z) coordinates
    # to the (i,j,k) axes.

    # For example, the first row of P (which contains pbest in column ibest)
    # determines the way the i axis points relative to the anatomical
    # (x,y,z) axes.  If ibest is 2, then the i axis is along the y axis,
    # which is direction P2A (if pbest > 0) or A2P (if pbest < 0).

    # So, using ibest and pbest, we can assign the output code for
    # the i axis.  Mutatis mutandis for the j and k axes, of course.

    return (ibest*pbest, jbest*qbest, kbest*rbest)
end

@inline function _nval_other_zero(n, val)
    if n === 1
        return val, 0, 0
    elseif n === 2
        return 0, val, 0
    else
        return 0, 0, val
    end
end

