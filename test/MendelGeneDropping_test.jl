using MendelBase
using DataFrames
using MendelGeneDropping

srand(123)

@testset "founder_source" begin
    # Every slice is a different locus, with columns being different people and columns 
    # being different alleles.
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

    #1st locus (Rh) not xlinked, so each person (columns) should have 2 source
    @test bush_matrix[:, 1, 1] == [1, 2]
    @test bush_matrix[:, 2, 1] == [3, 4]
    @test bush_matrix[:, 3, 1] == [5, 6]
    @test bush_matrix[:, 4, 1] == [0, 0] #only 3 founders in bush family, so last 3 column are 0
    @test bush_matrix[:, 5, 1] == [0, 0]
    @test bush_matrix[:, 6, 1] == [0, 0]

    # 2nd locus (ABO) not xlinked, so source number is also 1 ~ 6.
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

    #1st locus (Rh) not xlinked, so each person (columns) should have 2 source
    @test clinton_matrix[:, 1, 1] == [1, 2]
    @test clinton_matrix[:, 2, 1] == [3, 4]
    @test clinton_matrix[:, 3, 1] == [0, 0] #only 2 are founders

    # 2nd locus (ABO) is not xlinked, so source number is also 1 ~ 6
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
    x = [0.7, 0.2, 0.1]
    y = [0.5, 0.5, 0.5]
    z = [NaN, Inf]
    z2 = ["hi there"]

    cat_1, cat_2, cat_3 = 0, 0, 0
    cat_4, cat_5, cat_6 = 0, 0, 0
    cat_7, cat_8, cat_9 = 0, 0, 0
    for i in 1:100000
        test = MendelGeneDropping.random_genotype(x, true, true)
        if test == [1 1] cat_1 += 1 end
        if test == [1 2] cat_2 += 1 end
        if test == [1 3] cat_3 += 1 end
        if test == [2 1] cat_4 += 1 end
        if test == [2 2] cat_5 += 1 end
        if test == [2 3] cat_6 += 1 end
        if test == [3 1] cat_7 += 1 end
        if test == [3 2] cat_8 += 1 end
        if test == [3 3] cat_9 += 1 end                
    end

    #since male && xlinked, must have [i i] with probability around 0.7, 0.2, 0.1
    @test round(cat_1/100000, 1) == 0.7   
    @test round(cat_2/100000, 1) == 0.0
    @test round(cat_3/100000, 1) == 0.0
    @test round(cat_4/100000, 1) == 0.0    
    @test round(cat_5/100000, 1) == 0.2
    @test round(cat_6/100000, 1) == 0.0
    @test round(cat_7/100000, 1) == 0.0    
    @test round(cat_8/100000, 1) == 0.0
    @test round(cat_9/100000, 1) == 0.1

    cat_1, cat_2, cat_3 = 0, 0, 0
    cat_4, cat_5, cat_6 = 0, 0, 0
    cat_7, cat_8, cat_9 = 0, 0, 0
    for i in 1:100000
        test = MendelGeneDropping.random_genotype(x, false, true)
        if test == [1 1] cat_1 += 1 end
        if test == [1 2] cat_2 += 1 end
        if test == [1 3] cat_3 += 1 end
        if test == [2 1] cat_4 += 1 end
        if test == [2 2] cat_5 += 1 end
        if test == [2 3] cat_6 += 1 end
        if test == [3 1] cat_7 += 1 end
        if test == [3 2] cat_8 += 1 end
        if test == [3 3] cat_9 += 1 end                
    end

    @test round(cat_1/100000, 2) == 0.49 #0.7 * 0.7
    @test round(cat_2/100000, 2) == 0.14 #0.7 * 0.2 ...etc
    @test round(cat_3/100000, 2) == 0.07
    @test round(cat_4/100000, 2) == 0.14  
    @test round(cat_5/100000, 2) == 0.04
    @test round(cat_6/100000, 2) == 0.02
    @test round(cat_7/100000, 2) == 0.07
    @test round(cat_8/100000, 2) == 0.02
    @test round(cat_9/100000, 2) == 0.01

    @test_throws(ArgumentError, MendelGeneDropping.random_genotype(y, false, false))
    @test_throws(ArgumentError, MendelGeneDropping.random_genotype(z, false, false))
    @test_throws(MethodError, MendelGeneDropping.random_genotype(z2, false, false))
end

@testset "choose_genotype" begin
    x = [0.7, 0.2, 0.1]
    y = [0.5, 0.5, 0.5]
    z = [NaN, Inf]
    z2 = ["hi there"]

    genotype_set = Set{Tuple{Int64,Int64}}()
    @test_throws(ArgumentError, MendelGeneDropping.choose_genotype(y, 
        genotype_set, false, false))
    @test_throws(ArgumentError, MendelGeneDropping.choose_genotype(z, 
        genotype_set, false, false))
    @test_throws(MethodError, MendelGeneDropping.choose_genotype(z2, 
        genotype_set, false, false))

    cat_1, cat_2, cat_3 = 0, 0, 0
    cat_4, cat_5, cat_6 = 0, 0, 0
    cat_7, cat_8, cat_9 = 0, 0, 0
    for i in 1:100000
        test = MendelGeneDropping.choose_genotype(x, genotype_set, false, false)
        if test == [1 1] cat_1 += 1 end
        if test == [1 2] cat_2 += 1 end
        if test == [1 3] cat_3 += 1 end
        if test == [2 1] cat_4 += 1 end
        if test == [2 2] cat_5 += 1 end
        if test == [2 3] cat_6 += 1 end
        if test == [3 1] cat_7 += 1 end
        if test == [3 2] cat_8 += 1 end
        if test == [3 3] cat_9 += 1 end                
    end

    # everything is random because genotype_set is empty
    @test round(cat_1/100000, 2) == 0.49 #0.7 * 0.7
    @test round(cat_2/100000, 2) == 0.14 #0.7 * 0.2
    @test round(cat_3/100000, 2) == 0.07
    @test round(cat_4/100000, 2) == 0.14  
    @test round(cat_5/100000, 2) == 0.04
    @test round(cat_6/100000, 2) == 0.02
    @test round(cat_7/100000, 2) == 0.07
    @test round(cat_8/100000, 2) == 0.02
    @test round(cat_9/100000, 2) == 0.01

    x = [0.7, 0.1, 0.1, 0.1]
    genotype_set = Set{Tuple{Int64,Int64}}([(1, 2), (1, 1), (3, 1), (2, 3)])
    #first order the categories: (1, 1) (1, 2) (2, 3) (3, 1)
    #P(category 1) = 0.7 * 0.7
    #P(category 2) = 0.7 * 0.1
    #P(category 3) = 0.1 * 0.1
    #P(category 4) = 0.1 * 0.7
    #thus after normalizing, we have 0.765 0.109 0.0156 0.109 

    cat_1, cat_2, cat_3 = 0, 0, 0
    cat_4, cat_5, cat_6 = 0, 0, 0
    cat_7, cat_8, cat_9 = 0, 0, 0
    for i in 1:100000
        test = MendelGeneDropping.choose_genotype(x, genotype_set, false, false)
        if test == [1 1] cat_1 += 1 end
        if test == [1 2] cat_2 += 1 end
        if test == [1 3] cat_3 += 1 end
        if test == [2 1] cat_4 += 1 end
        if test == [2 2] cat_5 += 1 end
        if test == [2 3] cat_6 += 1 end
        if test == [3 1] cat_7 += 1 end
        if test == [3 2] cat_8 += 1 end
        if test == [3 3] cat_9 += 1 end                
    end

    @test round(cat_1/100000, 2) == 0.76 # (1 1)
    @test round(cat_2/100000, 2) == 0.11 # (1 2)
    @test round(cat_3/100000, 2) == 0.0
    @test round(cat_4/100000, 2) == 0.0 
    @test round(cat_5/100000, 2) == 0.0
    @test round(cat_6/100000, 2) == 0.02 # (2 3)
    @test round(cat_7/100000, 2) == 0.11 # (3 1)
    @test round(cat_8/100000, 2) == 0.0
    @test round(cat_9/100000, 2) == 0.0
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

    (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
        person, locus, keyword, 1) #1st pedigree, with 6 people and 3 generations

    @test size(sampled_genotype) == (2, 6, 4)
    @test eltype(sampled_genotype) == Int64
    @test all(sampled_genotype[:, :, :] .< 3) #test if value is 1 or 2
    @test all(sampled_genotype[:, :, :] .> 0) #test if value is 1 or 2

    # For foudners, check if sampled_genotype matrix and the founder matrix coincides
    bush_matrix = MendelGeneDropping.founder_source(pedigree, person, locus, 1, 
        keyword["gene_drop_output"])
    @test size(source) == (2, 6, 4)
    @test eltype(source) == Int64
    @test all(source[:, 1:3, :] .== bush_matrix[:, 1:3, :])

    #
    # First test for bush family
    # First two locus not xlinked
    #
    for j in 1:2
        allele_1_mom, allele_1_dad = 0, 0 #4th column 1st row 
        allele_2_mom, allele_2_dad = 0, 0 #4th column 2nd row
        allele_3_mom1, allele_3_dad1 = 0, 0 #5th column 1st row. needed 4 of these because that
        allele_3_mom2, allele_3_dad2 = 0, 0 #person could have either of the 4 alleles from the 1st generation
        allele_4_mom, allele_4_dad = 0, 0 #5th column 2nd row
        allele_5_mom, allele_5_dad = 0, 0 #6th column 1st row
        allele_6_mom, allele_6_dad = 0, 0 #6th column 2nd row
        for i in 1:100000
            (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
                person, locus, keyword, 1)
            #4th column
            if source[1, 4, j] == 3 allele_1_mom += 1 end
            if source[1, 4, j] == 4 allele_1_dad += 1 end
            if source[2, 4, j] == 1 allele_2_mom += 1 end
            if source[2, 4, j] == 2 allele_2_dad += 1 end

            #5th column
            if source[1, 5, j] == 1 allele_3_mom1 += 1 end
            if source[1, 5, j] == 2 allele_3_dad1 += 1 end
            if source[1, 5, j] == 3 allele_3_mom2 += 1 end
            if source[1, 5, j] == 4 allele_3_dad2 += 1 end
            if source[2, 5, j] == 5 allele_4_mom += 1 end
            if source[2, 5, j] == 6 allele_4_dad += 1 end
            
            #6th column
            if source[1, 6, j] == 3 allele_5_mom += 1 end
            if source[1, 6, j] == 4 allele_5_dad += 1 end
            if source[2, 6, j] == 1 allele_6_mom += 1 end
            if source[2, 6, j] == 2 allele_6_dad += 1 end        
        end
        @test round(allele_1_mom/100000, 2) == 0.5
        @test round(allele_1_dad/100000, 2) == 0.5
        @test round(allele_2_mom/100000, 2) == 0.5
        @test round(allele_2_dad/100000, 2) == 0.5
        @test round(allele_3_mom1/100000, 2) == 0.25
        @test round(allele_3_dad1/100000, 2) == 0.25  
        @test round(allele_3_mom2/100000, 2) == 0.25
        @test round(allele_3_dad2/100000, 2) == 0.25 
        @test round(allele_5_mom/100000, 2) == 0.5
        @test round(allele_5_dad/100000, 2) == 0.5
        @test round(allele_6_mom/100000, 2) == 0.5
        @test round(allele_6_dad/100000, 2) == 0.5       
    end

    #next 2 locus are xlinked
    for j in 3:4
        allele_1_mom, allele_1_dad = 0, 0 #4th column 1st row 
        allele_2_mom, allele_2_dad = 0, 0 #4th column 2nd row
        allele_3_mom1, allele_3_dad1 = 0, 0 #5th column 1st row. needed 4 of these because that
        allele_3_mom2, allele_3_dad2 = 0, 0 #person could have either of the 4 alleles from the 1st generation
        allele_4_mom, allele_4_dad = 0, 0 #5th column 2nd row
        allele_5_mom, allele_5_dad = 0, 0 #6th column 1st row
        allele_6_mom, allele_6_dad = 0, 0 #6th column 2nd row
        for i in 1:100000
            (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
                person, locus, keyword, 1)
            #4th column
            if source[1, 4, j] == 3 allele_1_mom += 1 end
            if source[1, 4, j] == 4 allele_1_dad += 1 end
            if source[2, 4, j] == 1 allele_2_mom += 1 end
            if source[2, 4, j] == 2 allele_2_dad += 1 end

            #5th column
            if source[1, 5, j] == 1 allele_3_mom1 += 1 end
            if source[1, 5, j] == 2 allele_3_dad1 += 1 end
            if source[1, 5, j] == 3 allele_3_mom2 += 1 end
            if source[1, 5, j] == 4 allele_3_dad2 += 1 end
            if source[2, 5, j] == 5 allele_4_mom += 1 end
            if source[2, 5, j] == 6 allele_4_dad += 1 end
            
            #6th column
            if source[1, 6, j] == 3 allele_5_mom += 1 end
            if source[1, 6, j] == 4 allele_5_dad += 1 end
            if source[2, 6, j] == 1 allele_6_mom += 1 end
            if source[2, 6, j] == 2 allele_6_dad += 1 end        
        end
        @test round(allele_1_mom/100000, 2) == 0.5
        @test round(allele_1_dad/100000, 2) == 0.5
        @test round(allele_2_mom/100000, 2) == 1.0
        @test round(allele_2_dad/100000, 2) == 0.0
        @test round(allele_3_mom1/100000, 2) == 0.5
        @test round(allele_3_dad1/100000, 2) == 0.0  
        @test round(allele_3_mom2/100000, 2) == 0.25
        @test round(allele_3_dad2/100000, 2) == 0.25
        @test round(allele_5_mom/100000, 2) == 0.5
        @test round(allele_5_dad/100000, 2) == 0.5
        @test round(allele_6_mom/100000, 2) == 1.0
        @test round(allele_6_dad/100000, 2) == 0.0       
    end

    #
    # Now test for Clinton family.
    # 1st two locus are not xlinked, so probability half for all 
    #
    for j in 1:2 
        allele_1_mom, allele_1_dad = 0, 0
        allele_2_mom, allele_2_dad = 0, 0
        for i in 1:100000
            (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
                person, locus, keyword, 2)
            if source[1, 3, j] == 3 allele_1_mom += 1 end
            if source[1, 3, j] == 4 allele_1_dad += 1 end
            if source[2, 3, j] == 1 allele_2_mom += 1 end
            if source[2, 3, j] == 2 allele_2_dad += 1 end
        end
        @test round(allele_1_mom/100000, 2) == 0.5
        @test round(allele_1_dad/100000, 2) == 0.5
        @test round(allele_2_mom/100000, 2) == 0.5
        @test round(allele_2_dad/100000, 2) == 0.5
    end

    #last two locus are xlinked, so since 1st person is male, his allele will surely be inherited 
    for j in 3:4 
        allele_1_mom, allele_1_dad = 0, 0
        allele_2_mom, allele_2_dad = 0, 0
        for i in 1:100000
            (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
                person, locus, keyword, 2)
            if source[1, 3, j] == 3 allele_1_mom += 1 end
            if source[1, 3, j] == 4 allele_1_dad += 1 end
            if source[2, 3, j] == 1 allele_2_mom += 1 end
            if source[2, 3, j] == 2 allele_2_dad += 1 end
        end
        @test round(allele_1_mom/100000, 2) == 0.5
        @test round(allele_1_dad/100000, 2) == 0.5
        @test round(allele_2_mom/100000, 2) == 1.0
        @test round(allele_2_dad/100000, 2) == 0.0
    end

    # Compute frequency for each combination of pedigree and locus. 
    # This can be done by running lines 245~283 in MendelGeneDropping.jl manually
    # freq2 = [0.61, 0.39, 0.0] #ped 1 locus 1. This is actually locus Rh, so we call it freq2
    # freq1 = [9/11, 2/11, 0.0, 0.0] #ped 1 locus 2. This is actually locus ABO 
    # freq3 = [0.67, 0.33, 0.0] #ped 1 locus 3. This is locus Xg
    # freq4 = [0.3, 0.7, 0.0] #ped 1 locus 4. This is locus XSNP

    # child1_allele1, child1_allele2 = 0, 0
    # child2_allele1, child2_allele2 = 0, 0
    # child3_allele1, child3_allele1 = 0, 0

    # for i in 1:100000
    #     (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
    #         person, locus, keyword, 1)
    #     if source[1, 4, 1] == 3 child1_allele1 += 1 end
    #     if source[1, 4, 1] == 4 child1_allele2 += 1 end
    #     if source[1, 4, 2] == 1 child1_allele2 += 1 end
    # end
end


@testset "convert_sampled_genotype" begin
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
    separator = "/"

    # first test clinton family, with 3 people only
    (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
                person, locus, keyword, 2)
    converted = MendelGeneDropping.convert_sampled_genotype(locus, 
        sampled_genotype, separator, "Unordered")
    @test size(converted) == (3, 4)
    @test eltype(converted) <: AbstractString

    # locus.name is ["Rh", "ABO", "Xg", "XSNP"] because the PEdigreeFrame does 
    # not have SNP as a header. So all information of SNP in locusframe is not included.
    @test converted[:, 1] == ["D/d", "D/D", "D/d"] #1 = D, 2 = d... etc
    @test converted[:, 2] == ["A/A", "A/A", "A/A"] 
    @test converted[:, 3] == ["+/+", "-/+", "+/+"]
    @test converted[:, 4] == ["2/2", "1/2", "2/2"]

    # now test bush family
    (sampled_genotype, source) = MendelGeneDropping.simulate_genotypes(pedigree, 
            person, locus, keyword, 1)
    converted = MendelGeneDropping.convert_sampled_genotype(locus, 
        sampled_genotype, separator, "Unordered")
    @test size(converted) == (6, 4)
    @test eltype(converted) <: AbstractString
    @test converted[:, 1] == ["D/D", "d/D", "D/d", "d/D", "d/D", "d/D"] #1 = D, 2 = d... etc
    @test converted[:, 2] == ["A/A", "A/A", "A/A", "A/A", "A/A", "A/A"] 
    @test converted[:, 3] == ["-/-", "+/+", "+/+", "+/-", "+/+", "+/-"]
    @test converted[:, 4] == ["1/1", "2/1", "2/2", "2/1", "2/2", "2/1"]
end



@testset "basics & wrapper functions" begin
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

    gene_drop = MendelGeneDropping.genedropping_option(pedigree, person, locus, 
        locus_frame, phenotype_frame, pedigree_frame, keyword)

    @test size(gene_drop) == (26, 17) #26 = 2 * 13, since 2 repetitions for each person
    @test eltype(gene_drop) == Any

    #testing if non-gene data is specified correctly
    @test gene_drop[1, 1] == "Bush1"
    @test gene_drop[1, 2] == "George"
    @test isna(gene_drop[1, 3])
    @test isna(gene_drop[1, 4])
    @test gene_drop[5, 5] == "female"
    @test gene_drop[6, 6] == 0

    #now test if ABO locus is being dropped according to probability specified in locus frame
    #can do so by looking at subset of population, like those who are 100% europeans
    AA_count = 0
    AB_count = 0
    BA_count = 0
    BB_count = 0
    for i = 1:20000
        gene_drop = MendelGeneDropping.genedropping_option(pedigree, person, locus, 
            locus_frame, phenotype_frame, pedigree_frame, keyword)
        if gene_drop[1, 12] == "A/A" AA_count+=1 end
        if gene_drop[1, 12] == "A/B" AB_count+=1 end
        if gene_drop[1, 12] == "B/A" BA_count+=1 end
        if gene_drop[1, 12] == "B/B" BB_count+=1 end
    end
    A_freq = 0.27 / (0.27 + 0.06) #these are specified in locusFrame. 
    B_freq = 0.06 / (0.27 + 0.06) #if person not 100% european, must adjust accordingly
    @test AA_count + AB_count + BA_count + BB_count == 20000 #no other blood type 
    @test round(AA_count/20000, 2) == round(A_freq^2, 2) 
    @test round(AB_count/20000, 2) == round(A_freq * B_freq, 2)
    @test round(BA_count/20000, 2) == round(A_freq * B_freq, 2) 
    @test round(BB_count/20000, 2) == round(B_freq^2, 2)

    final_test = GeneDropping("genedropping Control.txt")
    @test final_test == nothing #if everything ran smoothly, "nothing" should be returned.
end

#current coverage = (158, 201) ≈ 78.6%





