# PLANEJAMENTO DA OPERAÇÃO A CURTÍSSIMO PRAZO PARA DUAS USINAS EM CASCATA

Neste trabalho abordamos uma formulação para resolver o problema da alocação de unidades hidrelétricas, comumente chamado de Unit Commitment (UC).

O modelo trabalhado objetiva-se em minimizar a quantidade de água utilizada para atender as restrições do sistema, sendo aqui aplicadas para um horizonte diário com discretização horária. A modelagem deste problema é baseada em uma pesquisa que determina o status (ligado/desligado) e o nível de geração de cada unidade geradora (UG) para atender às restrições da usina e UGs, expandindo sua formulação para uma cascata. 

Foi implementado em Julia, sendo executado os experimentos numéricos com dados reais de operações de duas usinas e comparados com a operação real. O modelo se comportou adequadamente para aquilo que se propôs a fazer, sendo que os resultados se mostraram promissores com relação à economia da água nos reservatórios.
