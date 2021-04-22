using Distributions

mutable struct Lineage
  fitness::Float64
  size::Int64
  state::Vector{Float64}
  mod::Bool
end

function average_fitness(population::Dict{Array{Float64,1}, Lineage})
#calculate average fitness of the population. also outputs population size
  popN = Int64
  popN = 0
  popW = Float64
  popW = 0
  for line in values(population)
    popN+=line.size
    popW+=line.size * line.fitness
  end
  return popW/popN
end

function pop_size(population::Dict{Array{Float64,1}, Lineage})
  popN = Int64
  popN = 0
  for line in values(population)
    popN+=line.size
  end
  return popN
end

function assay_var_freq(population::Dict{Array{Float64,1}, Lineage})
#calculates  frequency of the sign-variable lineage
  popN = Int64
  popN = 0
  modN = Int64
  modN = 0
  for line in values(population)
    popN+=line.size
    if line.mod modN+=line.size end
  end
  modF = Float64
  modF = modN/popN
  return modF
end

function assay_var_count(population::Dict{Array{Float64,1}, Lineage})
#calculates the number of sign-variable carriers
  popN = Int64
  popN = 0
  modN = Int64
  modN = 0
  for line in values(population)
    popN+=line.size
    if line.mod modN+=line.size end
  end
  modF = Float64
  modF = modN/popN
  return modN
end


function model1(population::Dict{Array{Float64,1}, Lineage}, wcost::Float64, wben::Float64)
#between-lineage sign variability
  costly_start=false
  ben_start=false
  new_population = Dict{Array{Float64,1}, Lineage}()
  #initialize new population dictionary
  for (key, line) in population
    if line.mod
        pheno_dis = rand(Multinomial(line.size, [0.01, 0.99, 0.0]))
        #generates #s of beneficial, deleterious, and neutral carriers. When N=1 picks one of the three.
        if pheno_dis[3] > 0
            if key in keys(new_population)
                new_population[key].size+=pheno_dis[3]
            else new_population[key] = Lineage(line.fitness, pheno_dis[3], key, line.mod)
            end
        end

        if pheno_dis[1] > 0
            ben_start=true
            newkey = copy(line.state)
            newkey[2] =wben
            if newkey in keys(new_population) new_population[newkey].size+=pheno_dis[1]
            else new_population[newkey] = Lineage(wben, pheno_dis[1], newkey, line.mod)
            end
        end

        if pheno_dis[2] > 0
            costly_start=true
            newkey = copy(line.state)
            newkey[2] =wcost
            if newkey in keys(new_population) new_population[newkey].size+=pheno_dis[2]
            else new_population[newkey] = Lineage(wcost, pheno_dis[2], newkey, line.mod)
            end
        end

    else
        if key in keys(new_population) new_population[key].size+=line.size
        else new_population[key] = Lineage(line.fitness, line.size, line.state, line.mod)
        end
    end
  end
return new_population, costly_start, ben_start
end

function model2(generations::Int64, wcost::Float64, wben::Float64, t::Int64)
#within-lineage time-dependent sign variability model
  if generations < t
      w_inv=wcost
  else
      w_inv=wben
  end
  w_res=exp(1.0)
  return w_inv, w_res
end

function model3(freq::Float64, wcost::Float64, wben::Float64, f::Float64)
#within-lineage frequency-dependent sign variability model
  if freq < f
      w_inv=wcost
  else
      w_inv=wben
  end
  w_res=exp(1.0)
  return w_inv, w_res
end

function re_fit(population::Dict{Array{Float64,1}, Lineage}, w_inv::Float64, w_res::Float64)
    new_population = Dict{Array{Float64,1}, Lineage}()
    for (key, line) in population
        if line.mod
            newkey = copy(key)
            newkey[2] = w_inv
            new_population[newkey] = Lineage(w_inv, line.size, line.state, line.mod)
        else
            newkey = copy(key)
            newkey[2] = w_res
            new_population[newkey] = Lineage(w_res, line.size, line.state, line.mod)
        end
    end
    return new_population
end

function wright_fisher_reproduction(population::Dict{Array{Float64,1}, Lineage}, N::Int64, popw::Float64)
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[] #initialize lineage list
  for (key,line) in population
    push!(proby_list, line.size/N * line.fitness/popw)
    #representation of a lineage in the next generations depends on its frequency and relative fitness
    push!(lineage_list, line)
  end
  new_counts_list = rand(Multinomial(N, proby_list))
  new_population = Dict{Array{Float64,1}, Lineage}()
  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].mod)
    end
  end
  return new_population
end


function simulate(pop_N::Int64, model::Int64, wcost::Float64, wben::Float64, t::Int64, f::Float64)
  population = Dict{Array{Float64,1}, Lineage}()
  wt_state = [0.0, exp(1.0)]
  population[wt_state] = Lineage(exp(1.0), pop_N-1, wt_state, false)
  mod_state = [1.0, exp(1.0)]
  population[mod_state] = Lineage(exp(1.0), 1, mod_state, true)
  if model == 1 population, costly_start, ben_start = model1(population, wcost, wben) end
  modf = assay_var_freq(population)
  popw = average_fitness(population)
  generations = Int64
  generations = 0
  fgen = Int64
  fgen = 0
  # survivors = Int64
  # survivors = 0
  rfcrit=false
  fcr=0.0
  #initialize generations counter

  while 0.0<modf<1.0
    if model == 2
        w_inv, w_res = model2(generations, wcost, wben, t)
        population=re_fit(population,w_inv, w_res)
    elseif model == 3
        w_inv, w_res = model3(modf, wcost, wben, f)
        population=re_fit(population,w_inv, w_res)
    end
    popw = average_fitness(population)
    population = wright_fisher_reproduction(population, pop_N, popw)
    modf = assay_var_freq(population)
    popw = average_fitness(population)
    #if generations == 50 survivors =assay_var_count(population) end
    if modf < f
    else
        if rfcrit
        else
            rfcrit=true
            fgen=generations
            fcr=modf
        end
    end
    generations+=1
  end
  return modf ==1.0, generations, rfcrit, fgen, fcr, popw
end

job_id = ARGS[1]
for N in [10, 13, 17, 23, 32, 42, 57, 76, 102, 137, 183, 245, 327, 438, 586, 784, 1049, 1404, 1877, 2511]
  time = Float64[]
  # survs = Int64[]
  ftime= Int64[]
  f_cr=Float64[]
  println(N)
  #time_to_mut list records all times of successful mutator hitchhiking
  for run = 1:10000000
    output = simulate(N, 3, exp(0.99),exp(1.1), 50, 0.15)
    #simulation of model 3, records whether threshold frequency has been reached and the realized x* for every replicate
    if output[1]
        push!(time, output[2])
    else
        push!(time, 0)
    end

    if output[3]
        push!(ftime, output[4])
        push!(f_cr, output[5])

    else
        push!(ftime, 0)
        push!(f_cr, 0.0)

    end
    # push!(survs, output[4])
  end

outfile = open(string(job_id,"time.csv"), "a")
write(outfile, join(time, ","), "\n")
close(outfile)
outfile = open(string(job_id,"ftime.csv"), "a")
write(outfile, join(ftime, ","), "\n")
close(outfile)
outfile = open(string(job_id,"fcrt.csv"), "a")
write(outfile, join(f_cr, ","), "\n")
close(outfile)

end
