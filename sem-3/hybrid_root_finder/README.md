# EP1 de numérico

## Grupo: Anahí, Letícia, Isabella

---

## Como rodar

O arquivo `main.py` é diretamente executável em Linux, bastando fazer
`./main.py` para executar uma função teste.

O input para esse arquivo é no seguinte formato:

`./main.py [number] [errors] [interval] [constants]...`

- `number` indica o número da função desejada:
    0. Função teste
    1. Queda de um corpo sujeito ao atrito do ar
    2. Altura do fio
    3. Intercessão da curva da borboleta com a cardioide
    4. Intercessão da curva da borboleta com a terceira
    5. Intercessão da curva da borboleta com a quarta
    6. Intercessão da cardioide com a terceira curva
    7. Intercessão da cardioide com a quarta curva
    8. Intercessão da terceira curva com a quarta
- `errors` são os erros absoluto e relativo (escreva zero caso deseje utilizar
  o epsilon de ponto flutuante)
- `interval` indica o intervalo onde desejamos buscar raízes
- `constants` indica a escolha de constantes desejada pelo usuário

## Estrutura do código

### `main`

O ponto de entrada é o módulo `main`, onde temos o tratamento de input, seleção da função desejada e, então, chamamos a classe de aproximação.

### `functions`

As funções estão dispostas em arquivos específicos a depender de seus tipos
(i.e. `functions.physics` possui as funções com significado físico, e
`functions.polar_curves` possui as curvas polares). O módulo `functions`,
então, possui todas as funções num array `FUNCTIONS`.

É interessante notar também as meta funções definidas aqui, que devem guardar
meta-dados sobre as funções utilizadas no programa, como intervalo onde
desejamos buscar raízes (pré-selecionado), valores padrão para constantes e uma
breve descrição da função.

### `approximate`

Aqui está o coração do código onde, de fato, temos acesso ao método de Dekker e
também todas as funções de erro.

### `point`

Trata-se apenas de um arquivo que contém a definição de um ponto. É utilizado em outros módulos.
