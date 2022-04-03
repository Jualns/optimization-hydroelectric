
#mon = [577.58, 577.57, 577.57, 577.56, 577.56, 577.55, 577.55, 577.54, 577.54, 577.53, 577.53, 577.52, 577.52, 577.51, 577.51, 577.50, 577.50, 577.49, 577.49, 577.48, 577.48, 577.47, 577.47, 577.46]

#---------------------
#printl("A")

using CSV
using DataFrames
import XLSX


if ~occursin("Finardi", pwd())
    try
        cd("Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi")
    catch
        nothing
    end
end

#---
#Adicionar tempo da viagem da água #falta fazer
#Adicionar manutenção programada (paradas obrigatórias) #falta transformar em uma matriz binaria, testando.jl
#Modificar as entradas do modelo, colocar em matriz #falta fazer
#"Comparação" com a operação real (gráficos e médias de rendimentos, volume), sDemanda_Qx #feito, Afluencia.jl, melhorar eixos horas -> dias
#---


#vazao_real, potencia_real, Y = demanda_afluencia(dia = dia,mes,ano) # Demanda Espora para o dia 2
#L1 = potencia_real[1:end,4]
#L1 = demanda_afluencia(10,2019, mes_todo = 1) #Demanda Espora para o mês todo
#L1 = convert(Array{Float64,1}, Demanda_hora(02,10,2019; tipo = 1))# Demanda Espora antiga
xf = XLSX.readxlsx("Dados/hidr.xlsx") #abrindo a hidr
tamanho_final = string(xf["hidr"])[38:end-1] #tamanho inteiro [35:end-1]  da tabela, "A1:FX203"
sh = xf[string("hidr!A3:",tamanho_final)] #pego a tabela hidr do A3 até o tamanho_final

if ~isfile("Dados/entrada.csv")
    println("Arquivo \"entrada.csv\" faltante")
    cod_usinas = [71, 72, 73, 74, 76, 77]
    index_usinas = findall(x-> x in cod_usinas, sh[:,1]) #acha as linhas das usinas que queremos

    num_usinas = length(cod_usinas)
    nome_usina = sh[index_usinas, 2]

    aux_subsis = sh[index_usinas, 3] #subsistema por usina
    subsis = []
    for k = 1:num_usinas
        push!(subsis, split(aux_subsis[k]," - ")[1])
    end

    v_min = sh[index_usinas, 10]
    v_max = sh[index_usinas, 9]
    n_unidades = sh[index_usinas, 52]
    tipo_perda = sh[index_usinas, 46]
    perda_valor = sh[index_usinas, 47]
    w_min = sh[index_usinas,48] #vazão minima historica
    w_max = sh[index_usinas,54] # Vazão efetiva
    g_min = [10. 0. 0. 200. 150. 170.] #feito na mão pelo caderno de curvas colinas
    g_max = sh[index_usinas,53] #PotEf(1)

    jus_0 = sh[index_usinas, 147]
    jus_1 = sh[index_usinas, 148]
    jus_2 = sh[index_usinas, 149]
    jus_3 = sh[index_usinas, 150]
    jus_4 = sh[index_usinas, 151]

    mon_0 = sh[index_usinas, 15]
    mon_1 = sh[index_usinas, 16]
    mon_2 = sh[index_usinas, 17]
    mon_3 = sh[index_usinas, 18]
    mon_4 = sh[index_usinas, 19]


    min_var_h = []
    max_var_h = []
    for k = 1:num_usinas
        f(x) = jus_0[k] + jus_1[k]*x + jus_2[k]*x^2 + jus_3[k]*x^3 + jus_4[k]*x^4
        g(x) = mon_0[k] + mon_1[k]*x + mon_2[k]*x^2 + mon_3[k]*x^3 + mon_4[k]*x^4

        push!(min_var_h, g(sh[index_usinas[k],10]) - f(sh[index_usinas[k],54]*n_unidades[k])) #min - max
        push!(max_var_h, g(sh[index_usinas[k],9]) - f(0))#sh[index_usinas[k],48])) #max - min
    end

    df = DataFrame(cod_usinas = [71, 72, 73, 74, 76, 77],
    sub_sis = subsis,
    nome_usina = sh[index_usinas, 2],
    v_vert = sh[index_usinas, 13],
    v_min = sh[index_usinas, 10],
    v_max = sh[index_usinas, 9],
    n_unidades = sh[index_usinas, 52],
    tipo_perda = sh[index_usinas, 46],
    perda_valor = sh[index_usinas, 47],
    w_min = sh[index_usinas,48],
    w_max = sh[index_usinas,54],
    g_min = [10., 0., 0., 200., 150., 170.],
    g_max = sh[index_usinas,53],
    pot_total = g_max.*sh[index_usinas,52],
    jus_0 = sh[index_usinas, 147],
    jus_1 = sh[index_usinas, 148],
    jus_2 = sh[index_usinas, 149],
    jus_3 = sh[index_usinas, 150],
    jus_4 = sh[index_usinas, 151],
    mon_0 = sh[index_usinas, 15],
    mon_1 = sh[index_usinas, 16],
    mon_2 = sh[index_usinas, 17],
    mon_3 = sh[index_usinas, 18],
    mon_4 = sh[index_usinas, 19],
    min_h = min_var_h,
    max_h = max_var_h,
    rendi_nominal = sh[index_usinas, 37]*100)

    CSV.write("Dados/entrada.csv",df)
    println("Arquivo \"entrada.csv\" criado")
end

if ~isfile("Dados/potencia_instalada_hidr.csv")
    println("Arquivo \"potencia_instalada_hidr.csv\" faltante")
    subsistemas = unique(sh[2:end,3])
    subsis = []
    for k = subsistemas
        push!(subsis,split(k," - ")[1])
    end

    pot_max_instalada = []
    for k = subsistemas
        index_aux = findall(x -> x == k, sh[1:end,3])
        push!(pot_max_instalada, sum(sh[index_aux,52].*sh[index_aux,53] .+ sh[index_aux,56].*sh[index_aux,57] .+ sh[index_aux,60].*sh[index_aux,61] .+sh[index_aux,64].*sh[index_aux,65]))
    end

    df1 = DataFrame(sub_sis = subsis, #número do subsistema
    sub_sis_nome = subsistemas, #nome completo como está na hidr
    pot_instalada = pot_max_instalada)

    CSV.write("Dados/potencia_instalada_hidr.csv",df1)
    println("Arquivo \"potencia_instalada_hidr.csv\" criado")
end
