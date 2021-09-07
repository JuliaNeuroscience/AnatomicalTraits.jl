module AnatomicalTraits

using Metadata
using SpatioTemporalTraits
using Static

export is_sagittal, is_coronal, is_axial, is_radiologic, is_neurologic, is_anatomical

# TODO we need a simple method for mapping `src` to `dst`
# function permutation(src, dst) end

include("orientation.jl")

const SagittalName = Union{StaticSymbol{:left_to_right},StaticSymbol{:right_to_left},StaticSymbol{:sagittal}}

# TODO document
""" is_sagittal """
is_sagittal(x::Symbol) = is_sagittal(static(x))
is_sagittal(::SagittalName) = static(true)
is_sagittal(::Type{<:SagittalName}) = static(true)
is_sagittal(@nospecialize(x)) = static(false)

const CoronalName = Union{StaticSymbol{:anterior_to_posterior},StaticSymbol{:posterior_to_anterior},StaticSymbol{:coronal}}
# TODO document
""" is_coronal """
is_coronal(x::Symbol) = is_coronal(static(x))
is_coronal(::CoronalName) = static(true)
is_coronal(::Type{<:CoronalName}) = static(true)
is_coronal(@nospecialize(x)) = static(false)

const AxialName = Union{StaticSymbol{:superior_to_inferior},StaticSymbol{:inferior_to_superior},StaticSymbol{:axial}}
# TODO document
""" is_axial """
is_axial(x::Symbol) = is_axial(static(x))
is_axial(::AxialName) = static(true)
is_axial(::Type{<:AxialName}) = static(true)
is_axial(@nospecialize(x)) = static(false)

struct AnatomicalCoordinates <: AbstractCoordinateSystem end

@inline function orientation_name(::AnatomicalCoordinates, x::Int)
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
    is_radiologic(x) -> Bool

Returns `true` if `x` is consistent with radiologic orientation.
"""
is_radiologic(x) = is_radiologic(orientation(x))
@inline function is_radiologic(x::Tuple{Symbol,Symbol,Symbol})
    getfield(x, 1, false) === :left_to_right &&
    getfield(x, 2, false) === :anterior_to_posterior &&
    getfield(x, 3, false) === :superior_to_inferior
end
@inline function is_radiologic(x::Tuple{Symbol,Symbol})
    d1 = getfield(x, 1, false)
    d2 = getfield(x, 2, false)
    if d1 === :left_to_right
        # sagittal x coronal
        return d2 === :anterior_to_posterior
    else  # coronal x axial
        return d1 === :anterior_to_posterior && d2 === :superior_to_inferior
    end
end

"""
    is_neurologic(x) -> Bool

Returns `true` if `x` is consistent with neurologic orientation.
"""
is_neurologic(x) = is_neurologic(orientation(x))
@inline function is_neurologic(x::Tuple{Symbol,Symbol,Symbol})
    getfield(x, 1, false) === :right_to_left &&
    getfield(x, 2, false) === :anterior_to_posterior &&
    getfield(x, 3, false) === :superior_to_inferior
end
@inline function is_neurologic(x::Tuple{Symbol,Symbol})
    d1 = getfield(x, 1, false)
    d2 = getfield(x, 2, false)
    if d1 === :right_to_left
        # sagittal x coronal
        return d2 === :anterior_to_posterior
    else  # coronal x axial
        return d1 === :anterior_to_posterior && d2 === :superior_to_inferior
    end
end

"""
    is_anatomical(x) -> Bool

Returns `true` if `x` corresponds to anatomical data.
"""
@inline is_anatomical(x) = is_anatomical(orientation(x))
@inline function is_anatomical(x::Symbol)
    dynamic(is_sagittal(x)) || dynamic(is_coronal(x)) || dynamic(is_axial(x))
end
@inline function is_anatomical(x::Tuple{Symbol,Symbol})
    is_anatomical(getfield(x, 1, false)) && is_anatomical(getfield(x, 2, false))
end
@inline function is_anatomical(x::Tuple{Symbol,Symbol,Symbol})
    is_anatomical(getfield(x, 1, false)) &&
    is_anatomical(getfield(x, 2, false)) &&
    is_anatomical(getfield(x, 3, false))
end

# TODO document
""" AnatomicalRegion """
struct AnatomicalRegion{M}
    label::String
    metadata::M
end

Metadata.metadata(x::AnatomicalRegion) = getfield(x, :metadata)

label(x::AnatomicalRegion) = getfield(x, :label)

@inline abbreviation(x::AnatomicalRegion) = getmeta(label, x, :abbreviation)

@inline children(x::AnatomicalRegion) = getmeta(x, :children, ())

@inline hemisphere(x::AnatomicalRegion) = getmeta(x, :hemisphere, "")


end
