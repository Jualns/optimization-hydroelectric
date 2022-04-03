using Dates, XLSX, Plots, DataFrames
using CSV

if ~occursin("Finardi", pwd())
    try
        cd("Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi")
    catch
        nothing
    end
end

pasta_dados = "Dados/Dessem/"
entdados = string(pasta_dados,"entdados.dat")
#df = DataFrame(XLSX.readtable(string(,"/Dados/Dessem/entdat.xlsx"), "entdados"))
#names(df)
#println(df[!, 1])
#io = open(string(pasta_dados,"entdados.dat"), "r")
#a = read(io)
#close(io)

if ~isfile("Dados/Dessem/tempo_viagem.csv")
    df = CSV.File(entdados; header=false) |> DataFrame! #abrindo entdados.dat
    aux = df[!, 1] #pegando as informações
    pesquisando = ["&  TEMPO DE VIAGEM","& TERMICAS                             "] #os espaços importam
    posic = findall(x -> x in pesquisando, aux) #achando a posição da tabela "Tempo de Viagem"
    posic = posic .+ [2,-2] #excluindo extremidades textuais
    tempo_viagem = aux[posic[1]:posic[2]] #pegando tabela inteira
    colunas = Symbol.(split(tempo_viagem[1])[2:end]) #pegando o nome das colunas
    tempo_viagem = split.(tempo_viagem[3:end]) #pegando só os dados e transformando em matrix
    tempo_viagem = convert.(Array{Any}, tempo_viagem) #transformando o tipo para Any
    df_temp = Array(DataFrame(tempo_viagem)) #transformando em DataFrame e depois em Array para transport
    array_temp = permutedims(df_temp) #transpondo
    df = DataFrame(array_temp) #transformando em DataFrame
    df = df[!,2:end] #retirando a primeira coluna (TVIAG)
    rename!(df, colunas) #renomeando as colunas
    CSV.write(string(pasta_dados,"tempo_viagem.csv"), df)
    println("Arquivo \"tempo_viagem.csv\" criado")
end

if ~isfile("Dados/Dessem/demanda_sub_h.csv")
    df = CSV.File(entdados; header=false) |> DataFrame! #abrindo entdados.dat
    aux = df[!, 1] #pegando as informações
    pesquisando = ["&   CARGA","&   DEMANDAS/CARGAS ESPECIAIS"] #os espaços importam
    posic = findall(x -> x in pesquisando, aux) #achando a posição da tabela "Tempo de Viagem"
    posic = posic .+ [2,-3] #excluindo extremidades textuais
    demanda_sub_h = aux[posic[1]:posic[2]] #pegando tabela inteira
    colunas = Symbol.(split(demanda_sub_h[1])[2:end]) #pegando o nome das colunas
    colunas[4], colunas[7] = :mi, :mf
    demanda_sub_h = split.(demanda_sub_h[3:end]) #pegando só os dados e transformando em matrix
    demanda_sub_h = convert.(Array{Any}, demanda_sub_h) #transformando o tipo para Any
    df_temp = Array(DataFrame(demanda_sub_h)) #transformando em DataFrame e depois em Array para transport

    if size(df_temp)[1]<9 #9 é o número de colunas, contando a coluna de "DP"
        df_temp = vcat(df_temp[1:6,1:end],zeros(2,size(df_temp)[2]),df_temp[end:end,1:end])# adicionando duas colunas faltantes
    end

    array_temp = permutedims(df_temp) #transpondo
    df = DataFrame(array_temp) #transformando em DataFrame
    df = df[!,2:end] #retirando a primeira coluna (TVIAG)
    rename!(df, colunas) #renomeando as colunas
    CSV.write(string(pasta_dados,"demanda_sub_h.csv"), df)
    println("Arquivo \"demanda_sub_h.csv\" criado")
end
