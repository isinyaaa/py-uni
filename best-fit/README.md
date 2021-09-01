# Desafio Data Machina
Desafio da entrevista do Data Machina

## Problema #1 #
Nesse problema tive que pesquisar sobre API's e garantir a implementação correta de um método que geraria o número desejado de fibonacci.

### Utilização
Pela minha experiência com API's, nesse caso seria suficiente prover o método `get(n)` que retorna um objeto `response`, com o <img src="https://render.githubusercontent.com/render/math?math=n-">ésimo número de fibonacci (na variável `result`) e o status do serviço na variável `status`, ou seja:

* `result`
Nos dá a resposta desejada para a API.
* `status`
Nos dá o status da API, no momento só pode ser "invalid input" ou "success!". 

> Não é mais necessário instanciar um objeto da classe `fibo` para utilizar os comandos descritos.

## Problema #2 #

### Itinerário
- [X] Abstrações para o problema
- [X] Importação de dados
- [X] Organização dos itens dados
- [X] Algoritmo de encaixe
- [X] Estimativas de valores
- [X] Teste simples

#### Abstrações
Comecei estruturando os dados com os quais vamos trabalhar:

* `vehicle`
Um veículo tem sua capacidade máxima e essa que resta, e para cada cálculo podemos utilizar seus métodos internos.

* `item`
Um item é qualquer coisa que se queira transportar, e também tem atributos de dimensões.

* `company`
Uma empresa de frete qualquer, que possui sua frota de veículos e faz cotações, dado um carro e um peso qualquer.

#### Importação de dados
Utilizei a biblioteca padrão do `pandas` para isso, o que nos dá a flexibilidade de trabalhar com planilhas comuns exportadas em `.csv`.

#### Organização dos itens dados
Para a organização pensei em alguma medida do "quão quadrado é um item", e para isso lembrei do problema de maximização de área, onde:
> Dado um objeto retangular, queremos maximizar sua área
Sendo a solução padrão do problema utilizando a derivada e notando que a área é maximizada para objetos quadrados, posso utilizar um `sorting` simples analisando a derivada da área em função do perímetro e de uma das dimensões.

#### Algoritmo de encaixe
Esse teste deve determinar qual seria o menor carro onde caberiam todos os itens desejados, e como seriam colocados os itens:
1. Aqui podemos eliminar veículos menores classificando objetos pela sua altura e peso, que são atributos mais simples.
1. Depois podemos encaixar cada objeto normalmente ou trocando *comprimento* por *largura* (i.e. girando o objeto).
1. De acordo com [esse estudo informal](https://www.david-colson.com/2020/03/10/exploring-rect-packing.html) podemos esperar que um algoritmo de *row packing* seja bastante eficiente.
1. Podemos então usar algo como *column packing*, onde primeiro fazemos algo equivalente ao *row packing* com altura e área para depois realizar um row packing normal.
   > Cada coluna de objetos pode ser representada pelas dimensões do maior objeto naquela coluna (imagine várias caixas empilhadas num caminhão, as dimensões da maior caixa determinam onde iniciaremos outra pilha).
1. Finalmente, utilizamos *backtracking* para testar todas as possibilidades razoáveis numa margem aceitável de perda de performance/recursos.

#### Estimativas de valores
Essa estimativa pode ser feita estimando a "possibilidade de lucro" de cada veículo que poderia acomodar a carga, além da massa dos objetos alocados ali dentro.
> Dessa forma podemos estimar o consumo de combustível, além de uma taxa operacional mínima.

### Utilização
A tarefa foi separado em dois arquivos e diversas classes, como dito acima, mas a utilização permanece simplificada:
1. Os dados podem ser importados através da função `importData(name)` onde `name` é o nome de um arquivo `.csv` qualquer.
   > Essa função devolve uma lista de `item`'s.
1. Cada veículo pode ser criado utilizando sua própria classe, que recebe as dimensões máximas que ele suporta na ordem (comprimento, altura, largura, peso).
1. Cada companhia recebe uma lista de veículos com os quais vai trabalhar.
1. A função `listPossibilities(comp, data)` recebe uma companhia e os itens para trabalhar, e devolve uma lista de possiveis viagens para cada carro, no formato:
   > `Lista[índice do carro i] = [{itens carregados}, peso total do frete]` onde cada carro é organizado por ordem de capacidade.
> Todas essas funções estão localizadas dentro de `deliveryTools.py`. O arquivo `delivey.py` junto com o `meow.csv` compõem um exemplo de utilização.
