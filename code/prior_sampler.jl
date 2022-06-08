#################### Sampling from the prior ####################

#################### Types and Constructors ####################




mutable struct PriorTune <: SamplerTune
  logf::Union{Function, Missing}
  scale::Union{Float64, Vector{Float64}}
  eligible::Vector{Symbol}
  proposal::SymDistributionType

  PriorTune() = new()

  function PriorTune(x::Vector, scale::Real, logf::Union{Function, Missing}, eligible::Vector{Symbol};
                   proposal::SymDistributionType=Normal)
    new(logf, Float64(scale), eligible, proposal)
  end

  function PriorTune(x::Vector, scale::Vector{T},
                  logf::Union{Function, Missing},
                  eligible::Vector{Symbol};
                  proposal::SymDistributionType=Normal) where {T<:Real}
    new(logf, convert(Vector{Float64}, scale), eligible, proposal)
  end
end

const PriorVariate = Sampler{PriorTune, T} where T

validate(v::PriorVariate{T}) where T<:AbstractArray{S} where S<:Union{Real, GeneralNode} = validate(v, v.tune.scale)

validate(v::PriorVariate{T}, scale::Float64) where T<:AbstractArray{S} where S<:Union{Real, GeneralNode} = v 

function validate(v::PriorVariate{T}, scale::Vector) where T<:AbstractArray{S} where S<:Union{Real, GeneralNode}
  n = length(v)
  length(scale) == n ||
    throw(ArgumentError("length(scale) differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################
"""
    Prior(params::ElementOrVector{Symbol},
                  scale::ElementOrVector{T}; args...) where {T<:Real})

Construct a `Sampler` object for Prior sampling. Parameters are assumed to be
continuous, but may be constrained or unconstrained.

Returns a `Sampler{PriorTune}` type object.

* `params`: stochastic node(s) to be updated with the sampler. Constrained parameters are mapped to unconstrained space according to transformations defined by the Stochastic `unlist()` function.

* `scale`: scaling value or vector of the same length as the combined elements of nodes `params` for the `proposal` distribution. Values are relative to the unconstrained parameter space, where candidate draws are generated.

* `args...`: additional keyword arguments to be passed to the `PriorVariate` constructor.
"""
function Prior(params::ElementOrVector{Symbol}; args...) where {T<:Real}
  tune = PriorTune(Float64[], 1.0, logpdf!, Symbol[])
  Sampler(params, tune, Symbol[], false)
end


#################### Sampling Function ####################


"""
    sample!(v::PriorVariate, logf::Function)

Draw one sample from the distribution associated with each parameter

Returns `v` updated with simulated values and associated tuning parameters.
"""
function sample!(v::PriorVariate{T}, logf::Function; kwargs...) where T<:AbstractArray{<:Real}
  x = [rand(logf.m.nodes[p].distr) for p in v.params]
  v[:] = x
  v
end
