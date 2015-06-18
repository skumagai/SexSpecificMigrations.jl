import ForwardPopGenSimulations
const fpgs = ForwardPopGenSimulations

using SexSpecificMigrations
const sm = SexSpecificMigrations
using Base.Test

core = CoreData()
for i = 1:10
    @test sm.nextstate!(core) == i
    @test core.state == i
end
sm.settmax!(core, 10)
@test core.tmax == 10
gdb = sm.db(core)
@test isa(gdb, fpgs.GeneDB)
@test sm.time(core) == 0
for (i, t) in enumerate(core)
    @test i == t
end
@test sm.time(core) == 10

chrs = LocusList(
    [(Autosome, 2),
     (Autosome, 3),
     (XChromosome, 1),
     (YChromosome, 2),
     (Mitochondrion, 1)]
)
@test [chr for chr in chrs] == chrs.val
chrs = LocusList(
    (Autosome, 2),
    (Autosome, 3),
    (XChromosome, 1),
    (YChromosome, 2),
    (Mitochondrion, 1)
)
@test length(chrs) == 5
@test [chr for chr in chrs] == chrs.val

popsize = PopulationSizes(1=>(Female=>4, Male=>5), 2=>(Male=>6, Female=>7))
mut = MutationRates(Female=>0.5, Male=>0.2)
rec = RecombinationRates(Male=>0.2, Female=>0.1)
propm = ProportionsOfMigrants(1=>(Female=>0.2, Male=>0.3), 2=>(Female=>0.01, Male=>0.02))
fitm = FitnessOfMigrants(1=>(Female=>0.9, Male=>1.0), 2=>(Female=>0.2, Male=>0.1))
fitt = FitnessOfForeignTraits(1=>(Female=>0.8, Male=>0.9), 2=>(Female=>0.7, Male=>0.6))
tt = TraitTransmissionRates(Mother=>(Daughter=>0.2, Son=>0.3), Father=>(Daughter=>0.4, Son=>0.5))
@test length(popsize) == 2
@test length(propm) == 2
@test length(fitm) == 2
@test length(fitt) == 2
@test popsize[1] == (4, 5)
@test popsize[2] == (7, 6)
@test mut[1] == 0.5
@test mut[2] == 0.2
@test rec[1] == 0.1
@test rec[2] == 0.2
@test propm[1] == (0.2, 0.3)
@test propm[2] == (0.01, 0.02)
@test fitm[1] == (0.9, 1.0)
@test fitm[2] == (0.2, 0.1)
@test fitt[1] == (0.8, 0.9)
@test fitt[2] == (0.7, 0.6)
@test tt[1] == (0.2, 0.3)
@test tt[2] == (0.4, 0.5)
traits = Array{Float64}(10000)
for i = eachindex(traits)
    traits[i] = sm.trait(tt[1][1], tt[2][1], 1, 2, 1)
end
@test 0.35 < countnz(traits .== 2) / 10000 < 0.45
for i = eachindex(traits)
    traits[i] = sm.trait(tt[1][1], tt[2][1], 2, 1, 1)
end
@test 0.15 < countnz(traits .== 2) / 10000 < 0.25
for i = eachindex(traits)
    traits[i] = sm.trait(tt[1][1], tt[2][1], 2, 2, 1)
end
@test 0.47 < countnz(traits .== 2) / 10000 < 0.57

llist   = LocusList((Autosome, 1), (XChromosome, 1), (YChromosome, 1), (Mitochondrion, 1))
popsize = PopulationSizes(1=>(Female=>2, Male=>2), 2=>(Female=>2, Male=>2))
mut     = MutationRates(Female=>0.1, Male=>0.2)
rec     = RecombinationRates(Female=>0.2, Male=>0.3)
propm   = ProportionsOfMigrants(1=>(Female=>0.3, Male=>0.2), 2=>(Female=>0.2, Male=>0.1))
fitm    = FitnessOfMigrants(1=>(Female=>0.8, Male=>0.9), 2=>(Female=>0.7, Male=>0.6))
fitt    = FitnessOfForeignTraits(1=>(Female=>0.5, Male=>0.6), 2=>(Female=>0.6, Male=>0.7))
tt      = TraitTransmissionRates(Mother=>(Daughter=>0.9, Son=>0.2), Father=>(Daughter=>0.3, Son=>0.9))
params = ModelParameters(llist, popsize, mut, rec, propm, fitm, fitt, tt)
pops = sm.createpopulations(params)
@test length(pops) == 2
sexes = [Female, Male]
muts = [0.1, 0.2]
recs = [0.5, 0.5, 0.5, 0.5]
chrts = [Autosome, XChromosome, YChromosome, Mitochondrion]
for deme in eachindex(pops)
    for i in 1:2
        # For some reason, I can't compare SexSpecificPopulation to another with equality,
        # so I resorted to compare those field-by-field.

        @test sm.get(pops[deme][i], :deme) == deme
        @test size(pops[deme][i].data) == (2, 2 * 4)
        @test length(pops[deme][i]) == 2
        @test sm.get(pops[deme][i], :number_of_loci) == 4
        @test sm.get(pops[deme][i], :sex) == sexes[i]
        @test sm.get(pops[deme][i], :mutation_rate) == muts[i]
        for j = 1:sm.get(pops[deme][i], :number_of_loci)
            @test sm.get(pops[deme][i], :recombination_rate, j) == recs[j]
            @test sm.get(pops[deme][i], :type_of_locus, j) == chrts[j]
        end
        @test sm.get(pops[deme][i], :proportion_of_migrant) == propm[deme][i]
        @test sm.get(pops[deme][i], :fitness_of_foreign_trait) == fitt[deme][i]
        @test sm.get(pops[deme][i], :fitness_of_migrant) == fitm[deme][i]
        @test sm.get(pops[deme][i], :transmission_rate_of_trait, Daughter) == tt[i][1]
        @test sm.get(pops[deme][i], :transmission_rate_of_trait, Son) == tt[i][2]
        for j in 1:length(pops[deme][i])
            @test pops[deme][i][Val{:trait}, j] == deme
            pops[deme][i][Val{:trait}, j] = deme + 1
            @test pops[deme][i][Val{:trait}, j] == deme + 1
        end
        for j in eachindex(pops[deme][i])
            @test pops[deme][i][j] == 0
            pops[deme][i][j] = i
            @test pops[deme][i][j] == i
        end
        for j = 1:2, k = 1:2
            pops[deme][i][j, k] = 10 * j + k
            @test pops[deme][i][j, k] == 10 * j + k
        end
    end
end

pops = sm.createpopulations(params)
core = sm.CoreData()
@test core.state == 0
sm.initialize!(core, pops)
ans = (
    ([1 3 0 5 6 8 0 10; 2 4 0 0 7 9 0 0],
     [11 13 0 0 15 17 0 0; 12 0 14 0 16 0 18 0]),
    ([19 21 0 23 24 26 0 28; 20 22 0 0 25 27 0 0],
     [29 31 0 0 33 35 0 0; 30 0 32 0 34 0 36 0])
 )
@test pops[1][1].data == ans[1][1]
@test pops[1][2].data == ans[1][2]
@test pops[2][1].data == ans[2][1]
@test pops[2][2].data == ans[2][2]
@test maximum(collect(keys(sm.db(core)))) == 36
@test core.state == 36

for i = 1:2, j = 1:2
    for k in eachindex(pops[i][j].data)
        pops[i][j].data[k] = pops[i][j].data[k] == 0 ? 0 : pops[i][j].data[k] + 100
    end
end
for key in sort(collect(keys(core.db)))
    core.db[key + 100] = core.db[key]
    if core.db[key].parent.id == key
        core.db[key + 100].parent = core.db[key + 100]
    else
        core.db[key + 100].parent = core.db[core.db[key].parent.id + 100]
    end
end
for key in filter(x -> x < 100, keys(core.db))
    delete!(core.db, key)
end
core = sm.reinitialize!(core, pops)
@test core.state == 36
@test pops[1][1].data == ans[1][1]
@test pops[1][2].data == ans[1][2]
@test pops[2][1].data == ans[2][1]
@test pops[2][2].data == ans[2][2]

chrvec = Array{Int}(10000)
for i in eachindex(chrvec)
    chrvec[i] = sm.selectchr(1, 0.0)
end
@test all(chrvec[:] .== 1) == true
for i in eachindex(chrvec)
    chrvec[i] = sm.selectchr(1, 1.0)
end
@test all(chrvec[:] .== 2) == true
for i in eachindex(chrvec)
    chrvec[i] = sm.selectchr(1, 0.5)
end
@test 0.45 < countnz(chrvec[:] .== 1) / 10000 < 0.55

params2 = ModelParameters(
    chrs,
    PopulationSizes(1=>(Female=>2, Male=>2), 2=>(Female=>2, Male=>2)),
    MutationRates(Female=>0.1, Male=>0.2),
    RecombinationRates(Female=>0.2, Male=>0.3),
    ProportionsOfMigrants(1=>(Female=>0.3, Male=>0.2), 2=>(Female=>0.2, Male=>0.1)),
    FitnessOfMigrants(1=>(Female=>0.8, Male=>0.9), 2=>(Female=>0.7, Male=>0.6)),
    FitnessOfForeignTraits(1=>(Female=>0.5, Male=>0.6), 2=>(Female=>0.6, Male=>0.7)),
    TraitTransmissionRates(Mother=>(Daughter=>0.9, Son=>0.2), Father=>(Daughter=>0.3, Son=>0.9))
)
chrs = LocusList(
    [(Autosome, 2),
     (Autosome, 3),
     (XChromosome, 1),
     (YChromosome, 2),
     (Mitochondrion, 1)]
)

core = CoreData()
core.state = 1
core.t = 1
root = fpgs.GeneRecord(1, 1)
root.id = 1
insert!(core.db, root)

sm.mutate!(core, 1, 0.0)
@test core.state == 1
@test haskey(core.db, 2) == true
@test core.db[2].state == 1
@test core.db[2].epoch == 1
@test core.db[2].event == fpgs.Transmission()
core.t = 2
sm.mutate!(core, 2, 1.0)
@test core.state == 2
@test core.db[3].state == 2
@test core.db[3].epoch == 2
@test core.db[3].event == fpgs.Mutation()
@test core.db[3].parent == core.db[2]

function init(params)
    core = CoreData()
    sm.settmax!(core, 10)
    parpops = sm.createpopulations(params)
    sm.initialize!(core, parpops)
    chpops  = sm.createpopulations(params)
    parpops, chpops, core #, cr
end

ffa1, ffa2 = [0, 0], [0, 0]
mfa1, mfa2 = [0, 0], [0, 0]
ffx1, ffx2 = [0, 0], [0, 0]
fma1, fma2 = [0, 0], [0, 0]
mma1, mma2 = [0, 0], [0, 0]
for i = 1:10000
    parpops, chpops, core = init(params)
    sm.reproduce!(core, parpops[1], (1, 1), chpops[1][1], 1)
    gdb = sm.db(core)
    pid = gdb[chpops[1][1][1,1]].parent.id
    if pid == parpops[1][1][1,1]
        ffa1[1] += 1
    elseif pid == parpops[1][1][2,1]
        ffa1[2] += 1
    end
    pid = gdb[chpops[1][1][2,1]].parent.id
    if pid == parpops[1][2][1,1]
        mfa1[1] += 1
    elseif pid == parpops[1][2][2,1]
        mfa1[2] += 1
    end
    pid = gdb[chpops[1][1][1,2]].parent.id
    if pid == parpops[1][1][1,2]
        ffx1[1] += 1
    elseif pid == parpops[1][1][2,2]
        ffx1[2] += 1
    end
    pid = gdb[chpops[1][1][2,2]].parent.id
    @test pid == parpops[1][2][1,2]
    pid = gdb[chpops[1][1][1,4]].parent.id
    @test pid == parpops[1][1][1,4]
    @test chpops[1][1][1,3] == chpops[1][1][2,3] == chpops[1][1][2,4] == 0

    sm.reproduce!(core, parpops[2], (2, 2), chpops[1][1], 2)
    pid = gdb[chpops[1][1][1,5]].parent.id
    if pid == parpops[2][1][1,5]
        ffa2[1] += 1
    elseif pid == parpops[2][1][2,5]
        ffa2[2] += 1
    end
    pid = gdb[chpops[1][1][2,5]].parent.id
    if pid == parpops[2][2][1,5]
        mfa2[1] += 1
    elseif pid == parpops[2][2][2,5]
        mfa2[2] += 1
    end
    pid = gdb[chpops[1][1][1,6]].parent.id
    if pid == parpops[2][1][1,6]
        ffx2[1] += 1
    elseif pid == parpops[2][1][2,6]
        ffx2[2] += 1
    end
    pid = gdb[chpops[1][1][2,6]].parent.id
    @test pid == parpops[2][2][1,6]
    pid = gdb[chpops[1][1][1,8]].parent.id
    @test pid == parpops[2][1][1,8]
    @test chpops[1][1][1,7] == chpops[1][1][2,7] == chpops[1][1][2,8] == 0

    sm.reproduce!(core, parpops[1], (1, 1), chpops[1][2], 1)
    pid = gdb[chpops[1][2][1,1]].parent.id
    if pid == parpops[1][1][1,1]
        fma1[1] += 1
    elseif pid == parpops[1][1][2,1]
        fma1[2] += 1
    end
    pid = gdb[chpops[1][2][2,1]].parent.id
    if pid == parpops[1][2][1,1]
        mma1[1] +=1
    elseif pid == parpops[1][2][2,1]
        mma1[2] += 1
    end
    pid = gdb[chpops[1][2][2,3]].parent.id
    @test pid == parpops[1][2][2,3]
    @test chpops[1][2][2,2] == chpops[1][2][1,3] == chpops[1][2][1,4] == chpops[1][2][2,4] == 0

    sm.reproduce!(core, parpops[2], (2, 2), chpops[1][2], 2)
    pid = gdb[chpops[1][2][1,5]].parent.id
    if pid == parpops[2][1][1,5]
        fma2[1] += 1
    elseif pid == parpops[2][1][2,5]
        fma2[2] += 1
    end
    pid = gdb[chpops[1][2][2,5]].parent.id
    if pid == parpops[2][2][1,5]
        mma2[1] += 1
    elseif pid == parpops[2][2][2,5]
        mma2[2] += 1
    end
    pid = gdb[chpops[1][2][2,7]].parent.id
    @test pid == parpops[2][2][2,7]
    @test chpops[1][2][2,6] == chpops[1][2][1,7] == chpops[1][2][1,8] == chpops[1][2][2,8] == 0
end

@test sum(ffa1) == sum(ffa2) == sum(mfa1) == sum(mfa2) == sum(ffx1) == sum(ffx2) == 10000
@test sum(fma1) == sum(fma2) == sum(mma1) == sum(mma2) == 10000
@test 0.45 <= ffa1[1] / sum(ffa1) <= 0.55
@test 0.45 <= mfa1[1] / sum(mfa1) <= 0.55
@test 0.45 <= ffx1[1] / sum(ffx1) <= 0.55
@test 0.45 <= ffa2[1] / sum(ffa2) <= 0.55
@test 0.45 <= mfa2[1] / sum(mfa2) <= 0.55
@test 0.45 <= ffx2[1] / sum(ffx2) <= 0.55
@test 0.45 <= fma1[1] / sum(mfa1) <= 0.55
@test 0.45 <= mma1[1] / sum(mma1) <= 0.55
@test 0.45 <= fma2[1] / sum(mfa2) <= 0.55
@test 0.45 <= mma2[1] / sum(mma2) <= 0.55

parpops, chpops, core = init(params)
sm.initialize!(core, pops)

data = Array{Int}(16)
sm.listgenes!(data, parpops, 1)
@test data == [1, 2, 6, 7, 11, 12, 15, 16, 19, 20, 24, 25, 29, 30, 33, 34]
data = Array{Int}(12)
sm.listgenes!(data, parpops, 2)
@test data == [3, 4, 8, 9, 13, 17, 21, 22, 26, 27, 31, 35]
data = Array{Int}(4)
sm.listgenes!(data, parpops, 3)
@test data == [14, 18, 32, 36]
sm.listgenes!(data, parpops, 4)
@test data == [5, 10, 23, 28]
