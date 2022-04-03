
#o if a seguir verifica se o Julia está na pasta correta para buscar os arquivos necessários
if ~occursin("Finardi", pwd())
    try
        cd("Dropbox/LAC/Alocação de unidades/Curtissimo Cascata Finardi") #direciona o Julia para a pasta correta
    catch
        nothing
    end
end

#os if's a seguir verificam a existência dos arquivos necessários para rodar o modelo
if ~isfile("Dados/aflu_h.csv")
    include("aflu_h.jl")
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
if ~isfile("Dados/demanda_demanda.csv") || ~isfile("Dados/demanda_hidr.csv") #precisa ser o último dos if's
    include("entrada.jl")
end
#---------------------

using JuMP, Juniper, Ipopt, DataFrames
using CSV
using Dates
import XLSX
using Interpolations
using Plots
using LinearAlgebra

pyplot() #utilizando o pyplot

dia, mes, ano = 18,10,2019 #inserindo o dia para obtenção dos dados



#t1 = now()
optimizer = Juniper.Optimizer
nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
model = Model(optimizer_with_attributes(optimizer, "nl_solver" => nl_solver))#,"mip_solver"=>mip_solver))

grau = "2" #grau para coeficientes da curva colina
coefs_cc = CSV.File(string("Dados/Curva Colina/coefs_",grau,"grau.csv")) |> Tables.matrix
Hidr = CSV.File(string("Dados/entrada.csv")) |> DataFrame!

DefAnt = CSV.File(string("Dados/Dessem/deflu_h.csv")) |> DataFrame! #defluência anterior
#Dem_sub = CSV.File(string("Dados/Dessem/demanda_sub_h.csv")) |> DataFrame! #demanda

dem_dem = CSV.File(string("Dados/demanda_demanda.csv")) |> DataFrame! #demanda tipo demanda
dem_hidr = CSV.File(string("Dados/demanda_hidr.csv")) |> DataFrame! #demanda tipo hidr

temp_viag =  CSV.File(string("Dados/Dessem/tempo_viagem.csv")) |> DataFrame! #tempo de viagem
Aflu = CSV.File(string("Dados/aflu_h.csv")) |> DataFrame!

cod_usinas = collect(Iterators.flatten((71:74, 76:77))) #códigos da ONS das usinas presentes neste código
usinas = length(cod_usinas)

DefAnt_usinas = DefAnt[findall(x -> x in cod_usinas, DefAnt[!,:Jus]), 1:end]
temp_viag_usinas = temp_viag[findall(x -> x in cod_usinas, temp_viag[!,:Jus]), 1:end]

# 71 - Pot, liquida            # Sta Clara PR
# 72 - cc Constante/Fio d'agua # Fundão
# 73 - N temos cc/Reservatorio # Jordão
# 74 - Pot, liquida            # Foz do Areia/G.B. Munhoz
# 76 - Vaz, liquida            # Segredo
# 77 - Pot, bruta              # S. Santiago

usinas_geradoras = collect(Iterators.flatten((1:4, 6:usinas))) #pega todas usinas que geram algo, retira as que são somente reservatório
#usi_mon_natural = [[-1] [71] [72] [-1] [74] [73 76]] #usinas a montante no rio

#usi_mon_natural = Array{Any}(ones(6))#cod. usinas a montante no rio
#usi_mon_natural[6] = [3 5]
#usi_mon_natural[1:5] = [[0] [1] [2] [0] [4]]

##-------
#VARIÁVEIS MANUAIS
temp = 3 #períodos que o modelo irá rodar
tipo_aflu = "Media Movel" # = "Media Movel Scaled" #variável que seleciona entre média móvel ou média móvel escalonada
tipo_demanda = "hidr" # = dem #define o timpo de demanda utilizada durante o código
##-------


Vaz_aflu = Aflu[Aflu[!,:x3] .== tipo_aflu, 1:end]

v_vert = Hidr[!,"v_vert"] #vol min necessário para verter
v_min = Hidr[!,"v_min"] #
v_max = Hidr[!,"v_max"] #
H = 60.0 * 60.0 #contante de transformação de vazão [m³/s] para volume [m³]
n_unidades = Hidr[!,"n_unidades"]

#L é um dataframe ou seja pra pegar a demanda horária de uma usina com código X (Inteiro) se usa L[!,string(X)] ou L[!,Symbol(X)]
if tipo_demanda == "hidr"
    L = dem_hidr
else
    L = dem_dem
end

tipo_perda = Hidr[!,"tipo_perda"] #tipo de perda para altura líquida
perda_valor = Hidr[!,"perda_valor"] #valor utilizado na conta da altura líquida

w_min = Hidr[!,"w_min"]#limites da vazão turbinada
w_max = Hidr[!,"w_max"]#limites da vazão turbinada
w_min[3] = 10 #jordão turbina 10 por todo momento para não secar o resto do rio
w_max[3] = 10 #jordão turbina 10 por todo momento para não secar o resto do rio
g_min = Hidr[!,"g_min"]#limites da geração
g_max = Hidr[!,"g_max"]#limites da geração

I_max = ones(6)*3 #CHUTE
T_up = ones(6)*6 #tempo mínimo da máquina ligada CHUTE
T_down = ones(6)*2 #tempo mínimo da máquina desligada CHUTE
jus_0 = Hidr[!,"jus_0"]
jus_1 = Hidr[!,"jus_1"]
jus_2 = Hidr[!,"jus_2"]
jus_3 = Hidr[!,"jus_3"]
jus_4 = Hidr[!,"jus_4"]

mon_0 = Hidr[!,"mon_0"]
mon_1 = Hidr[!,"mon_1"]
mon_2 = Hidr[!,"mon_2"]
mon_3 = Hidr[!,"mon_3"]
mon_4 = Hidr[!,"mon_4"]

#---------------------
#max_h = Hidr[!,"max_h"]
#min_h = Hidr[!,"min_h"]


rendimento = Hidr[!,"rendi_nominal"] #"RENDIMENTO NOMINAL'


periodos_aflu = size(Aflu)[2] - 3 #tiro as 2 colunas de texto
Y = zeros(usinas, periodos_aflu)
O = zeros(usinas)
indice_usinas_montante = Dict() #dicionário que entra com o código da usina e retorna o índice i/k dela
δ = zeros(usinas)
for i in 1:usinas
    cod_aux = cod_usinas[i]
    if cod_aux in temp_viag_usinas[!,:Jus]
        δ[i] = temp_viag_usinas[temp_viag_usinas[!,:Jus] .== cod_aux,:hora][1]
    end

    posic_usina = findall(x -> x == cod_aux,  Vaz_aflu[!,:x2])[1]
    Y[i, 1:end] = Array(Vaz_aflu[posic_usina, 4:end])
    O[i] = sum(DefAnt_usinas[DefAnt_usinas[!,:Jus] .== cod_aux, :defluencia]) #vazões defluêntes anteriores

    get!(indice_usinas_montante, cod_aux, findall(x -> x in DefAnt_usinas[findall(x -> x == cod_aux, DefAnt_usinas[!,:Jus]),:Mont], cod_usinas))
end


#@variable(model, 1 ≥ x[ativo = 1:6] ≥ 0)
#OBS: Todas variáveis que contem n_unidades[5]*, não existem na entrada 5, logo o código não rodaria
@variable(model, q[i =1:usinas, t = 1:temp] ≥ 0)#, start = w_max[1]*n_unidades[1])#vazão turbinada total usina
@variable(model, s[i =1:usinas, t = 1:temp] ≥ 0)#, start = 0)#vazão vertida
@variable(model, v_max[i] ≥ v[i = 1:usinas, t = 0:temp] ≥ v_min[i])#, start = v_max[1])#volume
@variable(model, g[i =1:usinas, t = 1:temp, j = 1:n_unidades[i]]≥ 0)#, start = g_max[i]) #geração por turbina
@variable(model, hb[i=1:usinas, t = 1:temp] ≥ 0) #altura bruta, no código dá pra usar hl já que ele faz a transformação para hl ou deidxa em hb dependendo da cc
@variable(model, 1 ≥ η[i=1:usinas,t = 1:temp, j = 1:n_unidades[i]]≥ 0) #rendimento
@variable(model, w_max[i] ≥ w[i =1:usinas, t = 1:temp, j = 1:n_unidades[i]]≥ 0) #,start = w_max[i]) #vazão da turbina
@variable(model, u[i =1:usinas, t = 0:temp, j = 1:n_unidades[i]], Bin) #estado on off da turbina
@variable(model, ç[i =1:usinas, t = 0:temp, j = 1:n_unidades[i]], Bin) #ação de ligar
@variable(model, o[i =1:usinas, t = 0:temp, j = 1:n_unidades[i]], Bin) #ação de desligar
@variable(model, def[i=1:usinas, t = 1:temp] ≥ 0) #deficit
@variable(model, vaz_t[i = 1:usinas, t = 1:temp], start = 0.) #vazão transferida, positiva se a usina recebeu água, negativa se a usina passou água

#sta clara e foz do areia previsão de afluência por todo período
#fundão previsão de 2h, I_aflu = [10, 20, q+s[sta clara], ...]
#jordão previsão de 1h, I_aflu = [15, q+s[fundão], ...]

@NLexpression(model, mon[i = 1:usinas,
                        t = 1:temp],
        mon_4[i]*(v[i,t])^4 + mon_3[i]*(v[i,t])^3 +mon_2[i]*(v[i,t])^2 + mon_1[i]*(v[i,t]) +mon_0[i])#se der estranho pode ser aqui JORDAO

@NLexpression(model, jus[i = 1:usinas,
                        t = 1:temp],
        jus_4[i]*(q[i,t])^4 + jus_3[i]*(q[i,t])^3 +jus_2[i]*(q[i,t])^2 + jus_1[i]*(q[i,t]) +jus_0[i]) #se der estranho pode ser aqui JORDAO

#Restrição diz exatamente
#se i ∈ [1,4,5] então hl = (hb - perda_valor) ou (1-perda_valor/100)*hb
#caso contrário (i ∈ [2,3,6]) então hl = hb
@expression(model, hl[i = 1:usinas,
                    t = 1:temp],
        i in collect(Iterators.flatten((1:1, 4:5))) ? (2 - tipo_perda[i])*(1 - (perda_valor[i]/100.0))*hb[i,t] + (tipo_perda[i] - 1)*(hb[i,t] - perda_valor[i]) : hb[i,t])

@expression(model, pot[i=collect(Iterators.flatten((1, 4, 6))),
                    t = 1:temp,
                    j = 1:n_unidades[i]], #as cc da usina 1,4 e 6 são em pot
        9.80665*1e-3*rendimento[i]*hl[i,t]*w[i,t,j]) #usa a perda 1 se tipo perda == 1, análogo para 2


#Expressão que retorna a vazão transferida de Jordão para Segredo
#Se positiva -> Jordão recebe água de Segredo
#Se negativa -> Jordão envia água para Segredo
#Explicação: O túnel que liga as duas usinas funciona através da física dos vasos comunicantes
#@NLexpression(model, jord_segr[t = 1:temp], #vazão transferida de Jorão para Segredo
#        if mon[3,t] - mon[5,t] > 0
#            -72*sqrt(mon[3,t] - mon[5,t])
#        else
#            72*sqrt(-mon[3,t] + mon[5,t])
#        end
#        )
#Hidr[[3,6],[:max_h, :min_h, :v_min, :v_max]]
#collect(Iterators.flatten((1:4, 6:usinas)))
@expression(model,cc[i = 1:usinas,
                    t = 1:temp,
                    j = 1:n_unidades[i]],#[1,P,P*P,P*A,A,A*A]
    coefs_cc[i,1]+
    i == 5  ?  coefs_cc[i,2]*w[i,t,j] : coefs_cc[i,2]*g[i,t,j]  + #se i = 5 então w[] senão g[]
    i == 5  ?  coefs_cc[i,3]*w[i,t,j]^2 : coefs_cc[i,3]*g[i,t,j]^2 +
    i == 5  ?  coefs_cc[i,4]*w[i,t,j]*hl[t] : coefs_cc[i,4]*g[i,t,j]*hl[t] +
    coefs_cc[i,5]*hl[i,t]+
    coefs_cc[i,6]*hl[i,t]^2)

# Lembre que Segredo é com vazao_real
#@NLexpression(model,cc_segredo[i=5,t=1:temp, j=1:n_unidades[1]],#[1,P,P*P,P*A,A,A*A]
#    coefs_cc[i,1]+
#    coefs_cc[i,2]*w[i,t,j] +
#    coefs_cc[i,3]*w[i,t,j]^2 +
#    coefs_cc[i,4]*w[i,t,j]*hl[t] +
#    coefs_cc[i,5]*hl[i,t]+
#    coefs_cc[i,6]*hl[i,t]^2)

for k = 1:usinas
    for turbi = 1:n_unidades[k]
        fix.(ç[k,0,turbi],false , force = true)
        fix.(o[k,0,turbi],false , force = true)
        fix.(u[k,0,turbi],true , force = true)
    end
end

for t = 1:temp
    fix.(q[3,t],10,force = true)
end

for k = 1:usinas
    if ~(k in [3,6]) # para todas usinas menos Jordão (3) e Segredo (6)
        for t = 1:temp #em todos os tempos
            fix.(vaz_t[k,t], 0., force = true) #fixo vazão transferida = 0
        end
    end
    fix.(v[k,0], v_max[k] , force = true)## encontrar os dados iniciais de volume
end

@objective(model, Min, sum(q[i,t] + s[i,t] + 1e6*def[i,t] for i = 1:usinas, t = 1:temp))#REVER tem que multiplicar o deficit por alguma coisa!!

#se i = 5 => q[5,t] = 0.
#senão q[i,t] = sum(w[i,t,j] for j = 1:n_unidades[i])
@NLconstraint(model, turbinamento_usina[i = 1:usinas, t = 1:temp],
            q[i,t] == sum(w[i,t,j] for j = 1:n_unidades[i]))
            #i == 5 ? q[i,t] == 0. : q[i,t] == sum(w[i,t,j] for j = 1:n_unidades[i]))

#restrição do volume de vertimento(soleira), if v <= k, s = 0 elseif v > k -> s ∈ [0, s_max]
#@constraint(model, vol_soleira[i = 1:usinas, t = 1:temp], ########################################################
#            v[i,t] <= constante_por_usina ? s[i,t] = 0. : nothing)

#Se i = 3 então vaz_t[3,t] = vazão transferida de Jordão para Segredo
#Se i = 6 então vaz_t[6,t] = - vazão transferida de Jordão para Segredo
#Exemplo: Se Jordão manda 10 m3/s para Segredo então vaz_t[3,t] = -10 e vaz_t[6,t] = 10
vet_aux = zeros(usinas); vet_aux[3] = 1.; vet_aux[5] = -1. # usina 5 negativo pois se 3 recebe água então 5 perde e vice versa
@constraint(model,transferencia_vazao[i = 1:usinas, t = 1:temp],
            vaz_t[i,t] == vet_aux[i]*(mon[3,t] - mon[5,t])/abs(mon[3,t]- mon[5,t])*(-72*sqrt(abs(mon[3,t]- mon[5,t])))) #0 se i not in [3,5]
            #essa conta resulta na mesma operação que está abaixo:
            #if i == 3
            #    if mon[3,t] - mon[5,t] > 0
            #        vaz_t[i,t] == -72*sqrt(mon[3,t] - mon[5,t])
            #    else
            #        vaz_t[i,t] == 72*sqrt(-mon[3,t] + mon[5,t])
            #    end
            #elseif i == 5
            #    if mon[3,t] - mon[5,t] > 0
            #        vaz_t[i,t] == 72*sqrt(mon[3,t] - mon[5,t])
            #    else
            #        vaz_t[i,t] == -72*sqrt(-mon[3,t] + mon[5,t])
            #    end
            #else
            #    vaz_t[i,t] == 0.
            #end)

#usinas_geradoras = collect(Iterators.flatten((1:4, 6:usinas)) #pega todas usinas que geram algo, retira as que são somente reservatório
#Separando com if/elseif/else entre os rios
#Rio Jordão -> i = [1,2,3] observação: i = 6 também está no rio mas vamos considerar diferente
#Rio Iguaçu -> i = [4,5,6]
#Se i = 1 -> Estamos no rio Jordão
#Se i = 4 -> Estamos no rio Iguaçu
#Se i = 6 -> 1ᵃ usina que recebe os dois rios
@constraint(model,balanço_hidrico[i = 1:usinas, t = 1:temp],
            #v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]-Y[i,t]+vaz_t[i,t])) #o turbinado de jordão entra no incremental de Segredo e o vertido de jordão entra no incremental de S.Santiago
            if length(indice_usinas_montante[cod_usinas[i]]) == 0. #se é a primeira usina do rio
                v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]-Y[i,t]+vaz_t[i,t]))
            else
                if t > δ[i]
                    v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]- sum(q[i,t-Int(t-δ[i])] + s[i,t-Int(t-δ[i])] for i = indice_usinas_montante[cod_usinas[i]]) + vaz_t[i,t]))
                else
                    v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]- O[i] + vaz_t[i,t]))
                end
            end


            #if i == 1 (v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]-Y[i,t]+vaz_t[i,t]))
            #elseif i == 4 (v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]-Y[Int(i/2),t]+vaz_t[i,t]))
            #elseif i == 6 (v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]-q[i,t]-s[i,t]+vaz_t[i,t]))
            #else  (v[i,t] == v[i,t-1] - H.*(q[i,t]+s[i,t]-q[i,t]-s[i,t]+vaz_t[i,t]))
            #    sum(q[k,t] + s[k,t] for k = usi_mon_natural[i])

@constraint(model,balanço_hidrico_qx[ t = 1:temp],
             0 == (q[2,t]+s[2,t]-q[1,t]-s[1,t]))

@NLconstraint(model, atend_dem[i = collect(Iterators.flatten((1:2, 4:usinas))),t = 1:temp],
            sum(g[i,t,j] for j = 1:n_unidades[i]) == L[t,cod_usinas[i]])#rever o def[i,t] com o Tiago

@NLconstraint(model, geracao[i = collect(Iterators.flatten((1:4, 6:usinas))),t = 1:temp, j = 1:n_unidades[i]],
            9.80665*1e-3*η[i,t,j]*(1 -  perda_valor[i])*hb[t]*w[i,t,j] == g[i,t,j])

@NLconstraint(model, geracao_segredo[t = 1:temp, j = 1:n_unidades[5]],
            9.80665*1e-3*η[5,t,j]*hl*w[5,t,j]==g[5,t,j])

@NLconstraint(model, rendimento[i=collect(Iterators.flatten((1:4, 6:usinas))),t = 1:temp,j = 1:n_unidades[1]],
        cc[i,t,j]  == η[i,t,j])

@NLconstraint(model, rendimento_segredo[i=5,t = 1:temp,j = 1:n_unidades[1]],
        cc_segredo[i,t,j] == η[i,t,j])

@NLconstraint(model, altura_bruta[i=1:usinas,t = 1:temp],
            hb[i,t] - mon[i,t] + jus[i,t] == 0)

@constraint(model,lim_turbinamento_max[i = 1:usinas, t = 1:temp, j = 1:n_unidades[i]],
            w_max[i]*u[i,t,j] ≥ w[i,t,j])

@constraint(model,lim_turbinamento_min[i = 1:usinas,t = 1:temp, j = 1:n_unidades[i]],
            w[i,t,j] ≥ w_min[i]*u[i,t,j])

@constraint(model, lim_geracao_min[i = 1:usinas,t = 1:temp, j = 1:n_unidades[i]],
            g[i,t,j] ≥ g_min[i]*u[i,t,j])

@constraint(model, lim_geracao_max[i = 1:usinas, t = 1:temp, j = 1:n_unidades[i]],
            g_max[i]*u[i,t,j] ≥ g[i,t,j])

@constraint(model, start_up[i = 1:usinas,t = 1:temp, j = 1:n_unidades[i]],
            u[i,t,j]-u[i,t-1,j] == ç[i,t,j]-o[i,t,j])

@constraint(model, n_start_up_max[i = 1:usinas,j = 1:n_unidades[i]],
            sum(ç[i,t,j] for t = 1:temp) ≤ I_max[i] )

@constraint(model, n_star_up[i = 1:usinas,t = 1:temp,j = 1:n_unidades[i]],
            sum(ç[i,n,j] for n = max(1,t - T_up[i]):t ) ≤ u[i,t,j])

@constraint(model, n_shut_down[i = 1:usinas, j = 1:n_unidades[i], t = 1:temp],
            sum(o[i,n,j] for n = max(1,t - T_down[i]):t ) ≤ 1 - u[i,t,j])

optimize!(model)


error("AINDA NÃO")
if objective_value(model) > 0
    q1, s1, g1, w1 , v1, hb1 = value.(q), value.(s),value.(g),value.(w), value.(v), value.(hb)#, g, hb, w, n = value.(q), value.(s),value.(v), value.(g), value.(hb),value.(w), value.(η)
    #u1, I1, o1 = value.(ç),value.(I),value.(o)
    n = value.(η)
    println("----")
    println("Vazão turbinada usina $(q1)")
    println("----")
    println("Vazão vertida usina $(s1)")
    println("----")
    println("Volume usina $(v1)")
    println("----")
    println("Altura bruta $(hb1)")
    println("----")
    println("Geração:", g1)
    println("----")
    println("Vazão turbinada unidade:", w1)
    println("----")

    horas = 0:temp

    aux = []
    for k in 0:24
        push!(aux,v1[k])
    end

    volume = plot(horas, aux, xlabel = "Tempo (h)", ylabel ="Volume Útil (m³)" ,
        title = string("Volume Útil do Reservatório com Grau ",grau),  lw = 3, label = "Volume Útil")

    g_esp = ones(24,3)
    g_qx = ones(24,4)
    unidade = [3 4]
    for usi =1:2
        for t =1:24
            for uni =1:unidade[usi]
                if usi ==1
                    g_esp[t,uni] = g1[usi,t,uni]
                else
                    g_qx[t,uni] = g1[usi,t,uni]
                end
            end
        end
    end

    horas = 1:24
    geracao_esp = plot(horas, g_esp, xlabel = "Tempo (h)", ylabel ="Geração (MW/Médio)" ,
            title = string("Geração por Unidade na Usina 1 com Grau ",grau),label = ["Turbina 1" "Turbina 2" "Turbina 3"],  lw = 3)

    geracao_qx = plot(horas, g_qx, xlabel = "Tempo (h)", ylabel ="Geração (MW/Médio)",
            title = string("Geração por Unidade na Usina 2 com Grau ",grau),label = ["Turbina 1" "Turbina 2" "Turbina 3" "Turbina 4"],  lw = 3)

    w_esp = ones(24,3)
    w_qx = ones(24,4)
    n_esp = ones(24,3)
    n_qx = ones(24,4)
    for usi =1:2
        for t =1:24
            for uni = 1:unidade[usi]
                if usi ==1
                    w_esp[t,uni] = w1[usi,t,uni]
                    n_esp[t,uni] = n[usi,t,uni]
                else
                    w_qx[t,uni] = w1[usi,t,uni]
                    n_qx[t,uni] = n[usi,t,uni]
                end
            end
        end
    end

    vazao_turbi_esp = plot(horas, w_esp, xlabel = "Tempo (h)", ylabel ="Vazão Turbinada (m³/s)" ,
            title = string("Vazão turbinada por Unidade na Usina 1 com Grau ",grau),label = ["Turbina 1" "Turbina 2" "Turbina 3" ],  lw = 3)
    vazao_turbi_qx = plot(horas, w_qx, xlabel = "Tempo (h)", ylabel ="Vazão Turbinada (m³/s)" ,
            title = string("Vazão turbinada por Unidade na Usina 2 com Grau ",grau),label = ["Turbina 1" "Turbina 2" "Turbina 3" "Turbina 4"],  lw = 3)

    rendimento_esp = plot(horas, n_esp*100 , xlabel = "Tempo (h)", ylabel = "Rendimento (%)" ,
            title = string("Rendimento por Unidade na Usina 1 com Grau ",grau),label = ["Turbina 1" "Turbina 2" "Turbina 3" ],  lw = 3)

    #rend_qx = -0.18910719 .+ 0.11302677.*w_qx .- 0.0028301.*(w_qx.^2)
    rendimento_qx = plot(horas, n_qx*100 , xlabel = "Tempo (h)", ylabel = "Rendimento (%)" ,
                title = string("Rendimento por Unidade na Usina 2 com Grau ",grau),label = ["Turbina 1" "Turbina 2" "Turbina 3" "Turbina 4"],  lw = 3)

    d_esp = ones(24,2)
    d_qx = ones(24,2)
    for usi = 1:2
        for t = 1:24
            if usi ==1
                d_esp[t,1] = q1[usi,t]
                d_esp[t,2] = s1[usi,t]
            else
                d_qx[t,2] = s1[usi,t]
                d_qx[t,1] = q1[usi,t]
            end
        end
    end

    deflu_esp = plot(horas, [d_esp Y], xlabel = "Tempo (h)", ylabel ="Vazão (m³/s)" ,
            title = string("Defluência na Usina 1 com Grau ",grau),label = ["Turbinada" "Vertida" "Afluente"],  lw = 3)
    afluente = d_esp[:,1].+d_esp[:,2] #queixada
    deflu_qx = plot(horas, [d_qx afluente], xlabel = "Tempo (h)", ylabel ="Vazão (m³/s)" ,
            title = string("Defluência na Usina 2 com Grau ",grau),label = ["Turbinada" "Vertida" "Afluente"  ],  lw = 3)

    dem_esp = plot(horas, L1, xlabel = "Tempo (h)", ylabel ="Geração (MW/Médio)",
            title = string("Demanda na Usina 1 com Grau ",grau),label = "Demanda",  lw = 3)

    texto = string("cascata_espora")
    texto2 = string("cascata_queixada")
    pasta = string("Finardi/Todos")#/Cascata_sDemanda")
#----
    savefig(volume,string(pasta,"/vol_",texto,"_",grau,".svg"))
    savefig(geracao_esp,string(pasta,"/ger",texto,"_",grau,".svg"))
    savefig(vazao_turbi_esp,string(pasta,"/vaz_turbi_",texto,"_",grau,".svg"))
    savefig(deflu_esp,string(pasta,"/deflu_",texto,"_",grau,".svg"))
    savefig(rendimento_esp,string(pasta,"/rend_",texto,"_",grau,".svg"))
    savefig(dem_esp,string(pasta,"/dem_",texto,"_",grau,".svg"))

    savefig(geracao_qx,string(pasta,"/ger",texto2,"_",grau,".svg"))
    savefig(vazao_turbi_qx,string(pasta,"/vaz_turbi_",texto2,"_",grau,".svg"))
    savefig(deflu_qx,string(pasta,"/deflu_",texto2,"_",grau,".svg"))
    savefig(rendimento_qx,string(pasta,"/rend_",texto2,"_",grau,".svg"))
    println("Terminou")
end

df = DataFrame(horas = horas, Turb_esp = q1[1,1:end], Turb_qx = q1[2,1:end],
Vert_esp = s1[1,1:end], Vert_qx = s1[2,1:end],
Vol_oti = vol_oti[1:24,1], Vol_real = vol ,Hb_oti = hb1, Hb_real = h_dia,
g1_esp = g_esp[1:end,1], g2_esp = g_esp[1:end,2], g3_esp = g_esp[1:end,3],
g1_qx = g_qx[1:end,1], g2_qx = g_qx[1:end,2], g3_qx = g_qx[1:end,3], g4_qx = g_qx[1:end,4],
w1_esp = w_esp[1:end,1], w2_esp = w_esp[1:end,2], w3_esp = w_esp[1:end,3],
w1_qx = w_qx[1:end,1], w2_qx = w_qx[1:end,2], w3_qx = w_qx[1:end,3], w4_qx = w_qx[1:end,4],
n1_esp = n_esp[1:end,1], n2_esp = n_esp[1:end,2], n3_esp = n_esp[1:end,3],
n1_qx = n_qx[1:end,1], n2_qx = n_qx[1:end,2], n3_qx = n_qx[1:end,3], n4_qx = n_qx[1:end,4],
Dem_esp = L1,


CSV.write(string(pasta,"Dados_rodagem_",dia,"_",mes,"_",ano,".csv"), df)
t2 = now() - t1
#println("Tempo de execução do código: ",canonicalize(Dates.CompoundPeriod(t2)))
