using DataFrames
using CSV

#pasta_principal = "Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi"
if ~occursin("Finardi", pwd())
    try
        cd("Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi")

    catch
        nothing
    end
end

if ~isfile("Dados/entrada.csv") || ~isfile("Dados/potencia_instalada_hidr.csv")
    include("hidr.jl")
end
if ~isfile("Dados/Dessem/deflu_h.csv")
    include("deflu_h.jl")
end
if ~isfile("Dados/Dessem/demanda_sub_h.csv") || ~isfile("Dados/Dessem/tempo_viagem.csv")
    include("Entrada_dados.jl")
end



onde(codigo; cod_usinas = cod_usinas) = findall(x-> x == codigo, cod_usinas)[1]
#pasta_dados = "Dados/Dessem/"
cod_usinas = collect(Iterators.flatten((71:74, 76:77)))

#Defluência Anterior
DefAnt = CSV.File(string("Dados/Dessem/deflu_h.csv")) |> DataFrame!
index_def = findall(x -> x in cod_usinas, DefAnt[!,:Mont])
DefAnt = DefAnt[index_def,1:end]

#Pegando a Hidr (várias informações importantes das usinas estão aqui)
Hidr = CSV.File(string("Dados/entrada.csv")) |> DataFrame!
all_sub_sis = unique(Hidr[!,:sub_sis]) #Dict([(Hidr[k,:cod_usinas], Hidr[k,:sub_sis]) for k = 1:size(Hidr)[1]]) #cod_usina => subsistema
#gb = groupby(Hidr, :sub_sis)


Demanda_subsistema = Dict()
#Demanda por Subsistema
Dem_sub = CSV.File(string("Dados/Dessem/demanda_sub_h.csv")) |> DataFrame!
#demanda(sub_sistema) = Dem_sub[Dem_sub[!,:ss] .== sub_sistema,1:end]
#Dem_sub[Dem_sub[!,:ss] .== 2,1:end]
aux2 = Dem_sub[!,:di] .== 1 #dia inicial, chutei um dia qualquer

SubSis_usinas = Dict() #dicionário onde entra o subsistema e sai as usinas pertencentes a ele
for sub = all_sub_sis
    #get!(dicionario, variável, imagem), adiciona no dicionário variável => imagem
    #obs: não sobrescreve o valor caso variável já esteja no dicionário
    #get!(SubSis_usinas,sub,Array(gb[(sub_sis = 2,)][!,:cod_usinas]))
    get!(SubSis_usinas,sub, Hidr[Hidr[!,:sub_sis] .== sub, :cod_usinas])
    aux1 = Dem_sub[!,:ss] .== sub
    get!(Demanda_subsistema,sub, Dem_sub[aux1 .& aux2, :Demanda])
end

#Abor_horaria = false
#Demanda_usina = Dict()
#if Abor_horaria
#    for sub = all_sub_sis
#        usinas = SubSis_usinas[sub]
#        for usina = usinas
#            get!(Demanda_usina,usina
#        end
#    end
#end


pot_max_hidr = CSV.File(string("Dados/potencia_instalada_hidr.csv")) |> DataFrame!


max_hidr = Dict()
max_sub = Dict()
for sub = all_sub_sis
    get!(max_hidr, sub, pot_max_hidr[pot_max_hidr[!,:sub_sis] .== sub,:pot_instalada][1])
    get!(max_sub, sub, maximum(Demanda_subsistema[sub]))
end


#abordagem Hidr utiliza uma demanda menor por usina pois é > que o maximum da Demanda
if ~isfile("Dados/demanda_demanda.csv") #muda o limiter
    demanda_usina = Dict()
    #hidr = CSV.File(string("Dados/entrada.csv")) |> DataFrame!
    demanda_total  = 0.
    for sub = all_sub_sis
        contador = 0
        usinas_sub = SubSis_usinas[sub]
        limiter = max_sub[sub]
        for usina = usinas_sub
            pot_instalada_usina = Hidr[Hidr[!,:cod_usinas] .== usina,:pot_total][1]
            get!(demanda_usina, usina, Demanda_subsistema[sub] *pot_instalada_usina/limiter)
            #global nn += sum(demanda_usina[usina], dims = 2)
            if contador == 0
                global demanda_total = demanda_usina[usina]
            else
                global demanda_total += demanda_usina[usina]
            end
            contador += 1
        end
        get!(demanda_usina, -sub, demanda_total)
    end

    df_a = DataFrame(demanda_usina)
    CSV.write("Dados/demanda_demanda.csv", df_a)
    println("Arquivo \"demanda_demanda.csv\" criado")
end


if ~isfile("Dados/demanda_hidr.csv") #muda o limiter
    demanda_usina = Dict()
    #hidr = CSV.File(string("Dados/entrada.csv")) |> DataFrame!
    demanda_total  = 0.
    for sub = all_sub_sis
        contador = 0
        usinas_sub = SubSis_usinas[sub]
        limiter = max_hidr[sub]
        for usina = usinas_sub
            pot_instalada_usina = Hidr[Hidr[!,:cod_usinas] .== usina,:pot_total][1]
            get!(demanda_usina, usina, Demanda_subsistema[sub] *pot_instalada_usina/limiter)
            #global nn += sum(demanda_usina[usina], dims = 2)
            if contador == 0
                global demanda_total = demanda_usina[usina]
            else
                global demanda_total += demanda_usina[usina]
            end
            contador += 1
        end
        get!(demanda_usina, -sub, demanda_total)
    end

    df_b = DataFrame(demanda_usina)
    CSV.write("Dados/demanda_hidr.csv", df_b)
    println("Arquivo \"demanda_hidr.csv\" criado")
end
