using MendelBase
using DataFrames
using MendelGeneDropping

srand(123)

@testset "founder_source" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gene_drop_output"] = "Unordered"
    process_keywords!(keyword, "genedropping Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    bush_matrix = MendelGeneDropping.founder_source(pedigree, person, locus, 1, 
        keyword["gene_drop_output"])
    
    # NOTE: We have only 4 locus because the pedigree frame does not have "SNP" as header,
    #       So the SNP rows in the locus frame does not get read.  

    @test size(bush_matrix) == (2, 6, 4) #bush family have 6 people and 4 locus (instead of 5)
    @test eltype(bush_matrix) == Int64

    #1st locus (ABO) not xlinked, so each person (columns) should have 2 source
    @test bush_matrix[:, 1, 1] == [1, 2]
    @test bush_matrix[:, 2, 1] == [3, 4]
    @test bush_matrix[:, 3, 1] == [5, 6]
    @test bush_matrix[:, 4, 1] == [0, 0] #only 3 founders in bush family, so last 3 column are 0
    @test bush_matrix[:, 5, 1] == [0, 0]
    @test bush_matrix[:, 6, 1] == [0, 0]

    # 2nd locus (Rh) = 1st locus, since the 2nd locus is not xlinked, so source number should be the same.
    @test all(bush_matrix[:, :, 2] .== bush_matrix[:, :, 1])

    # 3rd locus (Xg) is xlinked, so males should only have 1 source while females have 2
    # NOTE: 2 is skipped, is that the expected behavior?
    @test bush_matrix[:, 1, 3] == [1, 1]
    @test bush_matrix[:, 2, 3] == [3, 4]
    @test bush_matrix[:, 3, 3] == [5, 5]
    @test all(bush_matrix[:, 4:6, 3] .== 0) 

    # 4th locus is XSNP which is xlinked (since SNP is skipped as explained in NOTE at the beginning)
    # so it behaves the same as the 3rd locus.
    @test all(bush_matrix[:, :, 4] .== bush_matrix[:, :, 3])


    # another test
    clinton_matrix = MendelGeneDropping.founder_source(pedigree, person, locus, 2, 
        keyword["gene_drop_output"])
    @test size(clinton_matrix) == (2, 3, 4) #clinton family have 3 people and 4 locus (instead of 5)
    @test eltype(clinton_matrix) == Int64

    #1st locus (ABO) not xlinked, so each person (columns) should have 2 source
    @test clinton_matrix[:, 1, 1] == [1, 2]
    @test clinton_matrix[:, 2, 1] == [3, 4]
    @test clinton_matrix[:, 3, 1] == [0, 0] #only 2 are founders

    # 2nd locus (Rh) = 1st locus, since the 2nd locus is not xlinked, so source number should be the same.
    @test all(clinton_matrix[:, :, 2] .== clinton_matrix[:, :, 1])
    
    # 3rd locus (Xg) is xlinked, so males should only have 1 source while females have 2
    # NOTE: 2 is skipped, is that the expected behavior?
    @test clinton_matrix[:, 1, 3] == [1, 1]
    @test clinton_matrix[:, 2, 3] == [3, 4]
    @test clinton_matrix[:, 3, 3] == [0, 0]

    # 4th locus is XSNP which is xlinked (since SNP is skipped as explained in NOTE at the beginning)
    # so it behaves the same as the 3rd locus.
    @test all(bush_matrix[:, :, 4] .== bush_matrix[:, :, 3])

    # need test cases for gene_drop_output == "Population"
end

@testset "random_genotype" begin
    x = [0.0, 0.0, 1.0]
    y = [0.7, 0.2, 0.1]
    z = [0.5, 0.5, 0.5]

    @test MendelGeneDropping.random_genotype(x, true, true) == [3 3]
    @test MendelGeneDropping.random_genotype(x, true, false) == [3 3]
    @test MendelGeneDropping.random_genotype(x, false, true) == [3 3]
    @test MendelGeneDropping.random_genotype(x, false, false) == [3 3]

    @test MendelGeneDropping.random_genotype(y, true, true) == [1 1]
    @test MendelGeneDropping.random_genotype(y, true, false) == [1 1]
    @test MendelGeneDropping.random_genotype(y, false, true) == [1 1]
    @test MendelGeneDropping.random_genotype(y, false, false) == [1 3]

    @test_throws(ArgumentError, MendelGeneDropping.random_genotype(z, false, false))
end

@testset "choose_genotype" begin
    
end

@testset "simulate_genotypes" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gene_drop_output"] = "Unordered"
    keyword["interleaved"] = true
    keyword["keep_founder_genotypes"] = false
    keyword["missing_data_pattern"] = "ExistingData" # Not yet implemented.
    keyword["missing_rate"] = 0.0
    keyword["repetitions"] = 1
    process_keywords!(keyword, "genedropping Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    (sampled_genotype2, source2) = MendelGeneDropping.simulate_genotypes(pedigree, 
        person, locus, keyword, 1) #only 3 pedigrees




end


@testset "convert_sampled_genotype" begin
    
end



@testset "basics & wrapper functions" begin
    
end