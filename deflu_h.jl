using Dates, XLSX, Plots, DataFrames
using CSV

#pasta_principal = "Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi"
if ~occursin("Finardi", pwd())
    try
        cd("Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi")

    catch
        nothing
    end
end

pasta_dados = "Dados/Dessem/"
deflant = string(pasta_dados,"deflant.dat")

df = CSV.File(deflant; header=false) |> DataFrame! #abrindo deflant.dat
##
aux = df[!, 1] #pegando as informações
posic = 4 #posição do começo da tabela
deflu_h = aux[posic:end] #pegando tabela inteira
colunas = Symbol.(split(deflu_h[1])[2:end]) #pegando o nome das colunas
colunas[6] = :mi
colunas[9] = :mf
deflu_h = split.(deflu_h[3:end]) #pegando só os dados e transformando em matrix
deflu_h = convert.(Array{Any}, deflu_h) #transformando o tipo para Any
df_temp = Array(DataFrame(deflu_h)) #transformando em DataFrame e depois em Array para transport
if size(df_temp)[1]<11 #11 é o número de colunas, contando a coluna de "DEFANT"
    df_temp = vcat(df_temp[1:8,1:end],zeros(2,size(df_temp)[2]),df_temp[end:end,1:end])# adicionando duas colunas faltantes
end
array_temp = permutedims(df_temp) #transpondo
df = DataFrame(array_temp) #transformando em DataFrame
df = df[!,2:end] #retirando a primeira coluna (DP)
rename!(df, colunas) #renomeando as colunas
CSV.write(string(pasta_dados,"deflu_h.csv"), df)
println("Arquivo \"deflu_h.csv\" criado")
