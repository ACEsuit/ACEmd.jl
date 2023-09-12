
@inline function ace_evaluate!(tmp, calc, R::AbstractVector, species::AbstractVector, species0 )
    return ACE1.evaluate!(tmp, calc, R, species, species0)
end

@inline function ace_evaluate!(B, tmp, calc, R::AbstractVector, species::AbstractVector, species0 )
    ACE1.evaluate!(B, tmp, calc, R, species, species0)
    return B
end

function ace_evaluate(calc, R::AbstractVector, species::AbstractVector, species0)
    tmp = ACE1.alloc_temp(calc, length(R))
    return ace_evaluate!(tmp, calc, R, species, species0)
end

function ace_evaluate(calc::ACE1.IPBasis, R::AbstractVector, species::AbstractVector, species0)
    B = ACE1.alloc_B(calc)
    tmp = ACE1.alloc_temp(calc, length(R))
    return ace_evaluate!(B, tmp, calc, R, species, species0)
end


@inline function ace_evaluate_d!(out, tmp, calc, R::AbstractVector, species::AbstractVector, species0)
    e = ACE1.evaluate_d!(out, tmp, calc, R, species, species0 )
    return e, tmp
end

function ace_evaluate_d(calc, R::AbstractVector, species::AbstractVector, species0)
    tmp = ACE1.alloc_temp_d(calc, length(R))
    return ace_evaluate_d!(tmp.dV, tmp, calc, R::AbstractVector, species::AbstractVector, species0)
end


function ace_evaluate_d(basis::ACE1.IPBasis, R::AbstractVector, species::AbstractVector, species0)
    dB  = ACE1.alloc_dB(basis, length(R))
    tmp = ACE1.alloc_temp_d(basis, length(R))
    ace_evaluate_d!(dB, tmp, basis, R, species, species0)
    return dB
end