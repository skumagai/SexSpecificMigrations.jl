using SexSpecificMigrations
using JLD

sttime = Dates.format(now(), "yyyy-m-d-H-M")
hostname = chomp(readall(`hostname`))
filename = "result.$hostname.$sttime.jld"

# One-locus per chromosomal type.
chrs = LocusList((Autosome, 1), (XChromosome, 1), (YChromosome, 1), (Mitochondrion, 1))
popsize = PopulationSizes(1=>(Female=>250, Male=>250), 2=>(Female=>250, Male=>250))
mutrates = MutationRates(Female=>0.000001, Male=>0.0000001)
recrates = RecombinationRates(Female=>0.1, Male=>0.1)
propm = ProportionsOfMigrants(
    1=>(Female=>0.0001, Male=>0.0001),
    2=>(Female=>0.0001, Male=>0.0001))
fitm = FitnessOfMigrants(1=>(Female=>0.8, Male=>0.9), 2=>(Female=>0.8, Male=>0.9))
fitt = FitnessOfForeignTraits(
    1=>(Female=>1.0, Male=>0.3),
    2=>(Female=>1.0, Male=>0.3))
trans = TraitTransmissionRates(
    Mother=>(Daughter=>0.9, Son=>0.6),
    Father=>(Daughter=>0.9, Son=>0.6))

params = ModelParameters(chrs, popsize, mutrates, recrates, propm, fitm, fitt, trans)

pops, genes, t = simulate(params, 5000, 30000, 3, 100)

@save filename pops genes t

