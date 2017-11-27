using MendelBase
using DataFrames
using MendelGeneDropping

srand(123)


@testset "founder_source" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gene_drop_output"] = "Unordered"
    # keyword["interleaved"] = true
    # keyword["keep_founder_genotypes"] = false
    # keyword["missing_data_pattern"] = "ExistingData" # Not yet implemented.
    # keyword["missing_rate"] = 0.0
    # keyword["repetitions"] = 1
    process_keywords!(keyword, "genedropping Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    matrix = MendelGeneDropping.founder_source(pedigree, person, locus, 1, 
        keyword["gene_drop_output"])
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

@testset "choose_genotype" begin
    
end

@testset "random_genotype" begin
    
end

@testset "convert_sampled_genotype" begin
    
end



@testset "basics & wrapper functions" begin
    
end