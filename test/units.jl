@testset "Unit conversion tests" begin
    folder = "./data/8-node-steady/"
    data = parse_data(folder)
    params, nominal_values = process_data!(data)
    make_english_units!(data, params, nominal_values)
    @test params[:is_english_units] == 1
    make_per_unit!(data, params, nominal_values)
    @test params[:is_per_unit] == 1
    make_si_units!(data, params, nominal_values)
    @test params[:is_si_units] == 1
end
