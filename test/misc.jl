
@testset "miscellaneous functions and corner cases" begin
    @test_throws DomainError S.sbml_export_sbo(
        A.Annotations("sbo" => ["123", "321"]),
    )

    @test S.sortunique([3, 2, 1, 2, 3, 2, 1, 2, 3]) == [1, 2, 3]
end
