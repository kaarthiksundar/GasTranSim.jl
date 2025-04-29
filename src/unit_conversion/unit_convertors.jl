function make_si_units!(
    data::Dict{String,Any},
    params::Dict{Symbol,Any},
    nominal_values::Dict{Symbol,Any},
)
    if params[:is_english_units] == 1
        english_to_si!(data, params, nominal_values)
        params[:is_si_units] = 1
        params[:is_english_units] = 0
        params[:is_per_unit] = 0
        return
    end
    if params[:is_per_unit] == 1
        pu_to_si!(data, params, nominal_values)
        params[:is_si_units] = 1
        params[:is_english_units] = 0
        params[:is_per_unit] = 0
        return
    end
    if params[:is_si_units] == 1
        return
    end
end

function make_per_unit!(
    data::Dict{String,Any},
    params::Dict{Symbol,Any},
    nominal_values::Dict{Symbol,Any},
)
    if params[:is_english_units] == 1
        english_to_si!(data, params, nominal_values)
        params[:is_si_units] = 1
        params[:is_english_units] = 0
        params[:is_per_unit] = 0
    end

    if params[:is_si_units] == 1
        si_to_pu!(data, params, nominal_values)
        params[:is_si_units] = 0
        params[:is_english_units] = 0
        params[:is_per_unit] = 1
        return
    end

    if params[:is_per_unit] == 1
        return
    end
end

function make_english_units!(
    data::Dict{String,Any},
    params::Dict{Symbol,Any},
    nominal_values::Dict{Symbol,Any},
)
    if params[:is_per_unit] == 1
        pu_to_si!(data, params, nominal_values)
        params[:is_si_units] = 1
        params[:is_english_units] = 0
        params[:is_per_unit] = 0
    end

    if params[:is_si_units] == 1
        si_to_english!(data, params, nominal_values)
        params[:is_si_units] = 0
        params[:is_english_units] = 1
        params[:is_per_unit] = 0
        return
    end

    if params[:is_english_units] == 1
        return
    end

end
