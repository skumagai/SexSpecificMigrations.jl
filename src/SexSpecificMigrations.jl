module SexSpecificMigrations

using ForwardPopGenSimulations

export Female,
       Male,
       Mother,
       Father,
       Daughter,
       Son,

       Autosome,
       XChromosome,
       YChromosome,
       Mitochondrion,

       CoreData,

       ModelParameters,
       LocusList,
       PopulationSizes,
       MutationRates,
       RecombinationRates,
       ProportionsOfMigrants,
       FitnessOfMigrants,
       FitnessOfForeignTraits,
       TraitTransmissionRates,

       simulate

type CoreData
    db::GeneDB
    state::Int
    t::Int
    tmax::Int
    CoreData() = new(GeneDB(), 0, 0)
end
nextstate!(d::CoreData) = d.state += 1
settmax!(d::CoreData, tmax) = d.tmax = tmax
db(d::CoreData) = d.db
time(d::CoreData) = d.t
Base.start(d::CoreData) = d.t = 0
Base.next(d::CoreData, state) = (d.t = state + 1; (d.t, d.t))
Base.done(d::CoreData, state) = state == d.tmax

immutable ChromosomeType
    val::Int
end
const Autosome = ChromosomeType(1)
const XChromosome = ChromosomeType(2)
const YChromosome = ChromosomeType(3)
const Mitochondrion = ChromosomeType(4)

immutable ForeignTrait end
immutable Migrant end

immutable SexType
    val::Int
end
const Female = SexType(1)
const Male = SexType(2)
# Mother, Father, Daughter, and Son are alias to Female and Male. These aliases are
# used in parameter specification to reduce confusion.
const Mother = SexType(1)
const Father = SexType(2)
const Daughter = SexType(1)
const Son = SexType(2)

immutable LocusList
    val::Array{Tuple{ChromosomeType, Int}, 1}
end
Base.call(::Type{LocusList}, val::Tuple{ChromosomeType, Int}...) = LocusList([val...])
Base.length(t::LocusList) = length(t.val)
Base.start(t::LocusList) = start(t.val)
Base.next(t::LocusList, s) = next(t.val, s)
Base.done(t::LocusList, s) = done(t.val, s)

abstract Parameter
abstract DemeAndSexSpecificParameter <: Parameter
abstract TwoGensSexSpecificParameter <: Parameter
abstract SexSpecificParameter <: Parameter
Base.getindex(p::Parameter, i::Integer) = p.data[i]

Base.length(i::DemeAndSexSpecificParameter) = length(i.data)

function Base.call{T<:DemeAndSexSpecificParameter, V<:Number}(
    ::Type{T}, vals::Pair{Int, NTuple{2, Pair{SexType, V}}}...)

    d = Dict([val[1] => Dict(val[2]) for val in vals])
    demes = sort(collect(keys(d)))
    T(tuple([(d[deme][Female], d[deme][Male]) for deme in demes]...))
end

function Base.call{T<:TwoGensSexSpecificParameter, V<:Number}(
    ::Type{T}, vals::Pair{SexType, NTuple{2, Pair{SexType, V}}}...)

    d = Dict([val[1] => Dict(val[2]) for val in vals])
    T(((d[Mother][Daughter], d[Mother][Son]), (d[Father][Daughter], d[Father][Son])))
end

function Base.call{T<:SexSpecificParameter, V<:Number}(
    ::Type{T}, vals::Pair{SexType, V}...)
    d = Dict(vals)
    T((d[Female], d[Male]))
end

immutable PopulationSizes <: DemeAndSexSpecificParameter
    # deme and sex-specific
    data::NTuple{2, NTuple{2, Int}}
    PopulationSizes(d::NTuple{2, NTuple{2, Int}}) = new(d)
end

immutable MutationRates <: SexSpecificParameter
    # sex-specific
    data::NTuple{2, Float64}
    MutationRates(d::NTuple{2, Float64}) = new(d)
end

immutable RecombinationRates <: SexSpecificParameter
    # sex-specific
    data::NTuple{2, Float64}
    RecombinationRates(d::NTuple{2, Float64}) = new(d)
end

immutable ProportionsOfMigrants <: DemeAndSexSpecificParameter
    # deme and sex-specific
    data::NTuple{2, NTuple{2, Float64}}
    ProportionsOfMigrants(d::NTuple{2, NTuple{2, Float64}}) = new(d)
end

abstract Fitness <: DemeAndSexSpecificParameter
immutable FitnessOfMigrants <: Fitness
    # deme and sex-specific
    data::NTuple{2, NTuple{2, Float64}}
    FitnessOfMigrants(d::NTuple{2, NTuple{2, Float64}}) = new(d)
end

immutable FitnessOfForeignTraits <: Fitness
    # deme and sex-specific
    data::NTuple{2, NTuple{2, Float64}}
    FitnessOfForeignTraits(d::NTuple{2, NTuple{2, Float64}}) = new(d)
end

immutable TraitTransmissionRates <: TwoGensSexSpecificParameter
    # sex-specific in parents and offspring
    data::NTuple{2, NTuple{2, Float64}}
    TraitTransmissionRates(d::NTuple{2, NTuple{2, Float64}}) = new(d)
end

function trait(momrate, poprate, momtrait, poptrait, deme)
    if deme == momtrait == poptrait
        deme
    elseif deme != momtrait && deme != poptrait
        rand() < 1.0 - (1.0 - momrate) * (1.0 - poprate) ? momtrait : deme
    elseif deme != momtrait
        rand() < momrate ? momtrait : deme
    else
        rand() < poprate ? poptrait : deme
    end
end

immutable ModelParameters
    chromosomes::LocusList
    population_sizes::PopulationSizes
    mutation_rates::MutationRates
    recombination_rates::RecombinationRates
    proportions_of_migrant::ProportionsOfMigrants
    fitness_of_migrants::FitnessOfMigrants
    fitness_of_foreign_traits::FitnessOfForeignTraits
    transmission_rates_of_trait::TraitTransmissionRates
end

immutable SexSpecificParameters
    deme::Int
    population_size::Int
    number_of_loci::Int
    sex::SexType
    mutation_rate::Float64
    recombination_rate::Array{Float64, 1}
    type_of_locus::Array{ChromosomeType, 1}
    proportion_of_migrant::Float64
    fitness_of_migrant::Float64
    fitness_of_foreign_trait::Float64
    transmission_rate_of_trait::NTuple{2, Float64}
end

function SexSpecificParameters(p::ModelParameters, deme::Integer, sex::SexType)
    sexidx = sex.val

    chrs = p.chromosomes
    recrates = vcat([[i == 1 ? 0.5 : recrate for i in 1:locus[2]] for locus in chrs]...)
    chrtypes = vcat([[locus[1] for _ in 1:locus[2]] for locus in chrs]...)

    length(recrates) == length(chrtypes) || error()

    nloci = length(chrtypes)

    SexSpecificParameters(
        deme,
        p.population_sizes[deme][sexidx],
        nloci,
        sex,
        p.mutation_rates[sexidx],
        recrates,
        chrtypes,
        p.proportions_of_migrant[deme][sexidx],
        p.fitness_of_migrants[deme][sexidx],
        p.fitness_of_foreign_traits[deme][sexidx],
        p.transmission_rates_of_trait[sexidx]
    )
end

type SexSpecificPopulation
    data::Array{Int, 2} # Chromosome is along row, and locus is along column.
    trait::Array{Int, 1}
    params::SexSpecificParameters
end

function SexSpecificPopulation(params::ModelParameters, deme::Int, sex::SexType)
    ps = SexSpecificParameters(params, deme, sex)
    nloci = length(ps.type_of_locus)
    popsize = ps.population_size
    SexSpecificPopulation(zeros(Int, 2, nloci * popsize), fill(deme, popsize), ps)
end

get(p::SexSpecificPopulation, s::Symbol) = getfield(p.params, s)
get(p::SexSpecificPopulation, s::Symbol, sex::SexType) = getfield(p.params, s)[sex.val]
get(p::SexSpecificPopulation, s::Symbol, i::Integer) = getfield(p.params, s)[i]
Base.setindex!(p::SexSpecificPopulation, val::Int, ::Type{Val{:trait}}, i::Integer) = p.trait[i] = val
Base.getindex(p::SexSpecificPopulation, ::Type{Val{:trait}}, i::Integer) = p.trait[i]
Base.length(p::SexSpecificPopulation) = p.params.population_size
Base.eachindex(p::SexSpecificPopulation) = eachindex(p.data)
Base.getindex(p::SexSpecificPopulation, i::Integer) = p.data[i]
Base.getindex(p::SexSpecificPopulation, c::Integer, l::Integer) = p.data[c, l]
Base.setindex!(p::SexSpecificPopulation, val::Int, i::Integer) = p.data[i] = val
Base.setindex!(p::SexSpecificPopulation, val::Int, c::Integer, l::Integer) = p.data[c, l] = val

type Population
    data::NTuple{2, SexSpecificPopulation}
end
Base.length(p::Population) = length(p.data)
Base.getindex(p::Population, i::Integer) = p.data[i]
Base.start(p::Population) = 0
Base.next(p::Population, state) = (state += 1; (p.data[state], state))
Base.done(p::Population, state) = length(p.data) == state
Base.eachindex(p::Population) = 1:length(p.data)

function createpopulations(params::ModelParameters)
    ndemes = length(params.population_sizes)
    pops = Array{Population}(ndemes)
    sexes = (Female, Male)
    for i = eachindex(pops)

        spops = Array{SexSpecificPopulation}(2)
        for (j, sex) = enumerate(sexes)
            spops[j] = SexSpecificPopulation(params, i, sex)
        end

        pops[i] = Population(tuple(spops...))
    end
    pops
end

function setinitialgenotype2!(core, p, i)
    for j = 1:2
        p[j,i] = insert!(db(core), GeneRecord(0, nextstate!(core)))
    end
    nothing
end

function setinitialgenotype1!(core, p, i, j)
    p[j,i] = insert!(db(core), GeneRecord(0, nextstate!(core)))
    nothing
end

function setinitialgenotype!(
    core::CoreData,
    p::SexSpecificPopulation,
    ind,
    locus,
    ::Type{Val{Female}},
    ::Union(Type{Val{Autosome}}, Type{Val{XChromosome}}))

    setinitialgenotype2!(core, p, toposition(p, ind, locus))
end

function setinitialgenotype!(
    core::CoreData,
    p::SexSpecificPopulation,
    ind,
    locus,
    ::Type{Val{Female}},
    ::Type{Val{YChromosome}})

    nothing
end

function setinitialgenotype!(
    core::CoreData,
    p::SexSpecificPopulation,
    ind,
    locus,
    ::Type{Val{Female}},
    ::Type{Val{Mitochondrion}})

    setinitialgenotype1!(core, p, toposition(p, ind, locus), 1)
end

function setinitialgenotype!(
    core::CoreData,
    p::SexSpecificPopulation,
    ind,
    locus,
    ::Type{Val{Male}},
    ::Type{Val{Autosome}})

    setinitialgenotype2!(core, p, toposition(p, ind, locus))
end

function setinitialgenotype!(
    core::CoreData,
    p::SexSpecificPopulation,
    ind,
    locus,
    ::Type{Val{Male}},
    ::Type{Val{XChromosome}})

    setinitialgenotype1!(core, p, toposition(p, ind, locus), 1)
end

function setinitialgenotype!(
    core::CoreData,
    p::SexSpecificPopulation,
    ind,
    locus,
    ::Type{Val{Male}},
    ::Type{Val{YChromosome}})

    setinitialgenotype1!(core, p, toposition(p, ind, locus), 2)
end

function setinitialgenotype!(
    core::CoreData,
    p::SexSpecificPopulation,
    ind,
    locus,
    ::Type{Val{Male}},
    ::Type{Val{Mitochondrion}})

    nothing
end

function initialize!(core, pops)
    for pop in pops
        for sexpop in pop
            nl = get(sexpop, :number_of_loci)
            sex = get(sexpop, :sex)
            popsize = length(sexpop)
            for j = 1:popsize, i in 1:nl
                chr = get(sexpop, :type_of_locus, i)
                setinitialgenotype!(core, sexpop, j, i, Val{sex}, Val{chr})
            end
        end
    end
end

function reinitialize!(oldcore, pops)
    oldgdb = db(oldcore)
    core = CoreData()
    gdb = db(core)
    # Preemptively map 0 to 0 as this value has a special meaning.
    smap = Dict{Int, Int}(0 => 0)
    smax = 1
    for pop in pops
        for sexpop in pop
            for idx in eachindex(sexpop.data)
                sexpop[idx] == 0 && continue
                state = oldgdb[sexpop[idx]].state
                if !haskey(smap, state)
                    smap[state] = smax
                    smax += 1
                end
                sexpop[idx] = insert!(gdb, GeneRecord(0, smap[state]))
            end
        end
    end
    core.state = maximum(keys(db(core)))
    core
end

function listgenes!(gids::AbstractArray, pops, locus, idx)
    for pop in pops
        for sexpop in pop
            chrtype = Val{get(sexpop, :type_of_locus, locus)}
            sex = get(sexpop, :sex)
            idx = listgenes!(gids, sexpop, locus, idx, chrtype, Val{sex})
        end
    end
    idx
end

listgenes!(gids::AbstractArray, pops, locus) = listgenes!(gids, pops, locus, 1)

listgenes!(gids::AbstractArray, sexpop, locus, idx, chrtype, sex) = idx

function listgenes!(gids::AbstractArray, sexpop, locus, idx, ::Type{Val{Autosome}}, sex)

    psize = length(sexpop)
    offset = get(sexpop, :number_of_loci)
    for i = 1:psize, chr = 1:2
        gids[idx] = sexpop[chr, (i - 1) * offset + locus]
        idx += 1
    end
    idx
end

function listgenes!(gids::AbstractArray, sexpop, locus, idx,
    ::Type{Val{XChromosome}}, ::Type{Val{Female}})

    psize = length(sexpop)
    offset = get(sexpop, :number_of_loci)
    for i = 1:psize, chr = 1:2
        gids[idx] = sexpop[chr, (i - 1) * offset + locus]
        idx += 1
    end
    idx
end

function listgenes!(gids::AbstractArray, sexpop, locus, idx,
    ::Type{Val{XChromosome}}, ::Type{Val{Male}})

    psize = length(sexpop)
    offset = get(sexpop, :number_of_loci)
    for i = 1:psize
        gids[idx] = sexpop[1, (i - 1) * offset + locus]
        idx += 1
    end
    idx
end


function listgenes!(gids::AbstractArray, sexpop, locus, idx,
    ::Type{Val{YChromosome}}, ::Type{Val{Male}})

    psize = length(sexpop)
    offset = get(sexpop, :number_of_loci)
    for i = 1:psize
        gids[idx] = sexpop[2, (i - 1) * offset + locus]
        idx += 1
    end
    idx
end

function listgenes!(gids::AbstractArray, sexpop, locus, idx,
    ::Type{Val{Mitochondrion}}, ::Type{Val{Female}})

    psize = length(sexpop)
    offset = get(sexpop, :number_of_loci)
    for i = 1:psize
        gids[idx] = sexpop[1, (i - 1) * offset + locus]
        idx += 1
    end
    idx
end

function reproduce!(core, parpop, ppos, chpop, cpos)
    # Maternally derived genes are always stored on the first row, whereas
    # paternally derived genes are on the second row.
    csex = get(chpop, :sex)
    for (sexpop, pos) in zip(parpop, ppos)
        psex = get(sexpop, :sex)
        nl = get(sexpop, :number_of_loci)
        cchr = psex.val
        chr = rand(1:2)
        for i = 1:nl
            chrtype = get(sexpop, :type_of_locus, i)
            geneid, chr = setgene!(core, Val{chrtype}, Val{psex}, Val{csex}, sexpop, pos, i, chr)
            chpop[cchr, toposition(chpop, cpos, i)] = geneid
        end
    end
end

selectchr(chr, rate) = rand() < rate ? 3 - chr : chr

function mutate!(core, id, mut)
    if rand() < mut
        transmit!(db(core), time(core), nextstate!(core), id)
    else
        transmit!(db(core), time(core), id)
    end
end

toposition(pop, ind, locus) = get(pop, :number_of_loci) * (ind - 1) + locus

setgene!(core, chr, psex, csex, pop, ind, locus, prevchr) = 0, prevchr

function setgene!(core, ::Type{Val{Autosome}}, psex, csex,
    pop, ind, locus, prevchr)

    pos = toposition(pop, ind, locus)
    chr = selectchr(prevchr, get(pop, :recombination_rate, locus))
    mutate!(core, pop[chr, pos], get(pop, :mutation_rate)), chr
end

function setgene!(core, ::Type{Val{XChromosome}}, ::Type{Val{Female}}, csex,
    pop, ind, locus, prevchr)

    pos = toposition(pop, ind, locus)
    chr = selectchr(prevchr, get(pop, :recombination_rate, locus))
    mutate!(core, pop[chr, pos], get(pop, :mutation_rate)), chr
end

function setgene!(core, ::Type{Val{Mitochondrion}}, ::Type{Val{Female}}, ::Type{Val{Female}},
    pop, ind, locus, prevchr)

    pos = toposition(pop, ind, locus)
    mutate!(core, pop[1, pos], get(pop, :mutation_rate)), 1
end

function setgene!(core, ::Type{Val{XChromosome}}, ::Type{Val{Male}}, ::Type{Val{Female}},
    pop, ind, locus, prevchr)

    pos = toposition(pop, ind, locus)
    mutate!(core, pop[1, pos], get(pop, :mutation_rate)), 1
end

function setgene!(core, ::Type{Val{YChromosome}}, ::Type{Val{Male}}, ::Type{Val{Male}},
    pop, ind, locus, prevchr)

    pos = toposition(pop, ind, locus)
    mutate!(core, pop[2, pos], get(pop, :mutation_rate)), 2
end

function evolve!(
    core::CoreData,
    parpops::Array{Population, 1},
    termon::Int,
    tclean::Int)

    ndemes = length(parpops) # basically fixed at 2
    nloci = get(parpops[1][1], :number_of_loci)
    loci = get(parpops[1][1], :type_of_locus)

    chpops = deepcopy(parpops)

    ncoals = [0 for _ = 1:nloci]
    lastcoals = [0 for _ = 1:nloci]

    # storage of gene ids. Autosome
    gidstore = Array{Array{Int, 1}, 1}(4)
    # Autosome
    gidstore[1] = zeros(
        Int,
        2 * sum((vcat([[length(pop) for pop in parpop] for parpop in parpops]...)))
    )
    # XChromosome
    gidstore[2] = zeros(
        Int,
        2 * sum([length(parpop[1]) for parpop in parpops]) +
        sum([length(parpop[2]) for parpop in parpops])
    )
    # YChromosome
    gidstore[3] = zeros(
        Int,
        sum([length(parpop[2]) for parpop in parpops])
    )
    # Mitochondrion
    gidstore[4] = zeros(
        Int,
        sum([length(parpop[1]) for parpop in parpops])
    )

    allgids = Array{Int}(sum(map(x -> length(gidstore[x.val]), loci)))

    # main loop of a simulation
    for gen = core
        for pop in chpops
            for sexpop in pop
                cdeme = get(sexpop, :deme)
                csex = get(sexpop, :sex)
                # i = toposition(sexpop)
                propmig = get(sexpop, :proportion_of_migrant)
                migfit = get(sexpop, :fitness_of_migrant)
                for org = 1:length(sexpop)
                    while true
                        # decide if an offspring is a native or migrant.
                        pdeme = cdeme
                        if rand() < propmig
                            rand() < migfit || continue
                            pdeme = 3 - cdeme
                        end

                        pp = parpops[pdeme]

                        # sample parents
                        mom = rand(1:length(pp[1]))
                        pop = rand(1:length(pp[2]))


                        # sample the trait of a offspring
                        momtrait = pp[1][Val{:trait}, mom]
                        poptrait = pp[2][Val{:trait}, pop]
                        ctrait = trait(
                            get(pp[1], :transmission_rate_of_trait, csex),
                            get(pp[2], :transmission_rate_of_trait, csex),
                            momtrait, poptrait, pdeme)

                        if pdeme != cdeme && ctrait != cdeme
                            tfit = get(sexpop, :fitness_of_foreign_trait)
                            rand() < tfit || continue
                        end

                        reproduce!(core, pp, (mom, pop), sexpop, org)

                        sexpop[Val{:trait}, org] = ctrait
                        break
                    end
                end
            end
        end

        parpops, chpops = chpops, parpops

        for (idx, locus) in enumerate(loci)
            lpos = locus.val
            listgenes!(gidstore[lpos], parpops, idx)
            anc = mrca(db(core), gidstore[lpos])
            if anc.epoch > lastcoals[idx]
                ncoals[idx] += 1
                lastcoals[idx] = gen
                info("turn over ", ncoals[lpos], " of locus ", locus, " on generation ", gen)
            end
        end
        if termon == minimum(ncoals)
            break
        end
        gen % tclean == 0 && cleandb!(core, parpops, allgids)
    end
    cleandb!(core, parpops, allgids)
    parpops, core

end

function cleandb!(core, pops, allgids)
    nl = get(pops[1][1], :number_of_loci)
    idx = 1
    for locus = 1:nl
        idx = listgenes!(allgids, pops, locus, idx)
    end
    clean!(db(core), allgids)
end

function simulate(params::ModelParameters, burnin::Int, t::Int, termon::Int, tclean::Int)
    core = CoreData()
    pops = createpopulations(params)

    initialize!(core, pops)
    settmax!(core, burnin)
    pops, core  = evolve!(core, pops, -1, tclean)

    core = reinitialize!(core, pops)
    settmax!(core, t)
    pops, core = evolve!(core, pops, termon, tclean)
    pops, db(core), time(core)
end

end # module
