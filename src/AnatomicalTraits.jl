module AnatomicalTraits

using Metadata
using SpatioTemporalTraits
using Static

include("orientation.jl")

const SagittalName = Union{StaticSymbol{:left_to_right},StaticSymbol{:right_to_left},StaticSymbol{:sagittal}}
is_sagittal(x::Symbol) = is_sagittal(static(x))
is_sagittal(::SagittalName) = static(true)
is_sagittal(::Type{<:SagittalName}) = static(true)
is_sagittal(@nospecialize(x)) = static(false)

const CoronalName = Union{StaticSymbol{:anterior_to_posterior},StaticSymbol{:posterior_to_anterior},StaticSymbol{:coronal}}
is_coronal(x::Symbol) = is_coronal(static(x))
is_coronal(::CoronalName) = static(true)
is_coronal(::Type{<:CoronalName}) = static(true)
is_coronal(@nospecialize(x)) = static(false)

const AxialName = Union{StaticSymbol{:superior_to_inferior},StaticSymbol{:inferior_to_superior},StaticSymbol{:axial}}
is_axial(x::Symbol) = is_axial(static(x))
is_axial(::AxialName) = static(true)
is_axial(::Type{<:AxialName}) = static(true)
is_axial(@nospecialize(x)) = static(false)


struct AnatomicalRegion{M}
    label::String
    metadata::M
end

Metadata.metadata(x::AnatomicalRegion) = getfield(x, :metadata)

label(x::AnatomicalRegion) = getfield(x, :label)

@inline abbreviation(x::AnatomicalRegion) = getmeta(label, x, :abbreviation)

@inline children(x::AnatomicalRegion) = getmeta(x, :children, ())

@inline hemisphere(x::AnatomicalRegion) = getmeta(x, :hemisphere, "")

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

end
