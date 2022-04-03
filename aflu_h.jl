#aflu_h.jl lê o arquivo Dados/Dados_hidro_ONS.csv, retirado do site da ONS com os dados de afluência
    #a discretização dessa vazão é diária para os últimos 5 anos, vou utilizar esse arquivo para obter algumas funções de vazão afluênte
        #a princípio penso em utilizar um modelo de média móvel
        #outra função seria pegando o max, min e média desses dados e construindo uma função
            #com o mesmo shape da demanda mas com esse min,max e média das vazões históricas
        #outra seria pegando o min, max e média e criando uma função randomica mas respeitando esses limites

using Dates, XLSX, Plots, DataFrames
using CSV
using Base.Sort

function n_menores(a, n, unico = 0) #pega os n menores valores do vetor a
    if unico == 0
        return sort(a; alg=Sort.PartialQuickSort(n))[1:n]
    end
    return sort(unique(a); alg=Sort.PartialQuickSort(n))[1:n]
end

function min_max(vec, min, max, matr = 1)
    v_min = minimum(vec, dims = 2)
    v_max = maximum(vec, dims = 2)
    if matr == 1
        return min .+ (max .- min).*(vec .- v_min)./(v_max .- v_min)
    end
    return min .+ (max - min)*(vec .- v_min)/(v_max - v_min)
end

#pasta_principal = "Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi"
if ~occursin("Finardi", pwd())
    try
        cd("Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi")

    catch
        nothing
    end
end

if ~isfile("Dados/ONS/Tratado_Afluencia.csv")
    error("Arquivo \"Tratado_Afluência.csv\" não encontrado")
end

#ma(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))] #média móvel padrão
ma_alterada(vs;n = 48) = [sum(@view vs[i:(i+(length(vs) - n))])/(length(vs) - n + 1) for i in 1:(length(vs)-(length(vs) - n))] #média móvel que acha a janela de dados que serão utilizados nas contas para retornar um vetor com N entradas

Vaz_aflu = CSV.File(string("Dados/ONS/Tratado_Afluencia.csv"), normalizenames=true) |> DataFrame!
quant_usinas = size(Vaz_aflu)[1]

array_vaz_aflu = Array(Vaz_aflu[!,3:end])
codigo_usinas = Vaz_aflu[!,2]

#2204 é o valor utilizado na média móvel como janela de dados utilizados para obter um vetor com 48 entradas no final
#Esse valor pode precisar de alterações dependendo do valor de entrada
all_vaz = ones(quant_usinas,48)
for usina = 1:quant_usinas
    all_vaz[usina,:] = ma_alterada(array_vaz_aflu[usina,1:end])#ma(Array(Vaz_aflu[usina,3:end]), 2204)
end

n = 5
menores = ones(quant_usinas)
maiores = ones(quant_usinas)
for k = 1:quant_usinas
    menores[k] = sum(n_menores(array_vaz_aflu[k,1:end],n,1))/n
    maiores[k] = sum(-n_menores(-array_vaz_aflu[k,1:end],n,1))/n
end

scaled_aflu = ones(size(all_vaz))
scaled_aflu = min_max(all_vaz, menores, maiores)
#aux = min_max.(array_vaz_aflu, menores, maiores)

nomes = []
reservatorio = []
for k = 1:2*quant_usinas
    if k <= quant_usinas
        push!(nomes, "Media Movel")
        push!(reservatorio, Vaz_aflu[k,:Reservatório])
    else
        push!(nomes, "Media Movel Scaled")
        push!(reservatorio, Vaz_aflu[k - quant_usinas,:Reservatório])
    end
end

df_vaz = DataFrame([reservatorio [codigo_usinas; codigo_usinas] nomes [all_vaz; scaled_aflu]])

CSV.write("Dados/aflu_h.csv", df_vaz)
