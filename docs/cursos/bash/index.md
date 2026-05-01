# Minicurso 1: Introdução ao Bash para a Bioinformática

!!! abstract "Resumo"
    Tutorial introdutório ao uso do terminal Bash aplicado à bioinformática. Cobre navegação no sistema de arquivos, visualização e manipulação de dados biológicos (logs, FASTA, Newick), uso de expressões regulares, loops e estruturas condicionais — ferramentas essenciais para quem trabalha com dados moleculares em servidores Linux.

## :material-target: Objetivos de aprendizagem

Ao final deste tutorial, você será capaz de:

- [ ] Navegar no sistema de arquivos de um servidor Linux
- [ ] Criar, mover, copiar e remover arquivos e diretórios
- [ ] Visualizar e inspecionar o conteúdo de arquivos biológicos
- [ ] Extrair informações com `grep`, `cut`, `sed` e pipes
- [ ] Usar expressões regulares para identificar padrões em dados
- [ ] Criar e utilizar variáveis, arrays e loops `for`/`while`
- [ ] Escrever estruturas condicionais `if` e `case`
- [ ] Automatizar tarefas repetitivas com arquivos FASTA e listas

## :material-clock-outline: Carga horária

**4 horas**

## :material-tools: Pré-requisitos

| Requisito | Nível esperado |
|-----------|---------------|
| Acesso a terminal Linux/macOS ou WSL | Básico |
| Conta em servidor remoto (SSH) | Desejável |
| Conhecimento prévio de bioinformática | Não necessário |

---

## Parte 1 — Primeiros Passos

### 1.1 Se encontrando no servidor

```bash
pwd
# Mostra o caminho completo do diretório atual
```

---

### 1.2 Listando arquivos e se localizando no servidor

```bash
ls          # Lista simples dos arquivos no diretório atual
ls *        # Lista todos os arquivos e exibe o conteúdo dos diretórios do primeiro nível
ls -lah     # Lista detalhada com arquivos ocultos e tamanhos legíveis

# Flags utilizadas:
# -l  →  Exibe detalhes em formato de lista (permissões, dono, tamanho, etc.)
# -a  →  Mostra arquivos ocultos (que começam com ".")
# -h  →  Mostra tamanhos em formato legível (KB, MB, GB)
```

| Flag | Significado |
|------|-------------|
| `-l` | Formato de lista detalhada |
| `-a` | Inclui arquivos ocultos |
| `-h` | Tamanhos legíveis por humanos |

---

### 1.3 Mostrar a hierarquia de pastas

```bash
tree  # Imprime a hierarquia de diretórios em formato de árvore
```

---

### 1.4 Mudança de diretórios

```bash
cd .    # Permanece no diretório atual
cd ..   # Sobe um nível na hierarquia
cd      # Vai para o diretório home (equivalente a cd ~)
cd ~    # Vai para o diretório home
cd /    # Vai para o diretório raíz
```

| Comando | Destino |
|---------|---------|
| `cd .` | Diretório atual |
| `cd ..` | Diretório pai |
| `cd ~` | Diretório home |
| `cd /` | Raíz do sistema |

---

### 1.5 Criação e manipulação de arquivos e pastas

```bash
nano arquivo0.txt   # Abre o editor de texto Nano
mkdir nova_pasta    # Cria um novo diretório
mv arquivo0.txt nova_pasta  # Move o arquivo para o diretório nova_pasta
```

---

### 1.6 Inspecionar o conteúdo de arquivos

```bash
wc -l nova_pasta/arquivo0.txt  # Conta o número de linhas
wc -w nova_pasta/arquivo0.txt  # Conta o número de palavras
wc -c nova_pasta/arquivo0.txt  # Conta o número de caracteres
```

| Flag | O que conta |
|------|-------------|
| `-l` | Linhas |
| `-w` | Palavras |
| `-c` | Caracteres (bytes) |

---

### 1.7 Copiar, renomear e remover arquivos

```bash
cp nova_pasta/arquivo0.txt backup.txt  # Copia o arquivo
mv backup.txt backup.renomeado         # Renomeia o arquivo
rm backup.renomeado                    # Remove o arquivo
```

### Checklist da Parte 1

- [ ] Executei `pwd` e entendi meu caminho no servidor
- [ ] Listei arquivos com `ls -lah` e interpretei as colunas
- [ ] Naveguei entre diretórios com `cd`
- [ ] Criei um arquivo com `nano` e um diretório com `mkdir`
- [ ] Copiei, renomeei e removi arquivos

---

## Parte 2 — Visualização de Dados

```bash
cd ~/1.visualizacao_de_dados
```

### 2.1 Visualização de conteúdo de arquivos

```bash
head -n 10 arquivo1.log   # Imprime as 10 primeiras linhas
tail -n 10 arquivo1.log   # Imprime as 10 últimas linhas

# Ver todo o conteúdo do arquivo:
cat  arquivo1.log   # Imprime todo o conteúdo de uma vez
less arquivo1.log   # Paginação interativa (use q para sair)
more arquivo1.log   # Paginação básica
```

| Comando | Comportamento |
|---------|---------------|
| `head -n N` | Primeiras N linhas |
| `tail -n N` | Últimas N linhas |
| `cat` | Todo o conteúdo |
| `less` | Paginação interativa |
| `more` | Paginação básica |

---

### 2.2 Pegando informações de interesse com pipes

```bash
grep -i "taxa" arquivo1.log | cut -d " " -f1,2 | sed "s/taxa/especies/g" > relatorio.info
```

**O que cada comando faz nessa pipeline:**

| Etapa | Comando | Função |
|-------|---------|--------|
| 1 | `grep -i "taxa" arquivo1.log` | Busca linhas com "taxa" (sem distinção de maiúsculas) |
| 2 | `cut -d " " -f1,2` | Seleciona as colunas 1 e 2 (separadas por espaço) |
| 3 | `sed "s/taxa/especies/g"` | Substitui "taxa" por "especies" em todo o resultado |
| 4 | `> relatorio.info` | Salva o resultado no arquivo (sobrescreve se existir) |

```bash
grep -i "Best-fit model" arquivo1.log >> relatorio.info
# >> adiciona ao final do arquivo sem apagar o conteúdo anterior
```

> **Dica:** `>` sobrescreve; `>>` acrescenta ao final.

---

### 2.3 Pegando espécies de um arquivo Newick (expressões regulares)

```bash
cd ~/1.visualizacao_de_dados/1.2_newicks
```

```bash
grep -E -o "[A-Z][a-z]*_[a-z]*:|[A-Z][a-z]*_[a-z]*_[a-z]*:" arquivo2.newick \
    | cut -d ":" -f1 > species.list
```

**Explicação dos padrões de expressão regular:**

| Padrão regex | Exemplo de match | Descrição |
|--------------|-----------------|-----------|
| `[A-Z][a-z]*_[a-z]*:` | `Homo_sapiens:` | Gênero + espécie + dois-pontos |
| `[A-Z][a-z]*_[a-z]*_[a-z]*:` | `Canis_lupus_familiaris:` | Gênero + espécie + subespécie + dois-pontos |

**Flags utilizadas:**

| Flag | Função |
|------|--------|
| `-E` | Ativa expressões regulares estendidas |
| `-o` | Imprime apenas o trecho que corresponde ao padrão |

---

### 2.4 Pegar o menor valor de bootstrap

```bash
grep -o "[0-9]*:" arquivo2.newick | cut -d ":" -f1 | sort -n | uniq
```

**Pipeline passo a passo:**

| Etapa | Comando | O que faz |
|-------|---------|-----------|
| 1 | `grep -o "[0-9]*:"` | Extrai números seguidos de `:` |
| 2 | `cut -d ":" -f1` | Remove o `:` final |
| 3 | `sort -n` | Ordena numericamente |
| 4 | `uniq` | Remove duplicatas consecutivas |

### Checklist da Parte 2

- [ ] Visualizei as primeiras e últimas linhas de um arquivo com `head` e `tail`
- [ ] Construí uma pipeline com `grep`, `cut` e `sed`
- [ ] Usei `>` para criar e `>>` para acrescentar a um arquivo
- [ ] Extrai nomes de espécies de um arquivo Newick com expressões regulares
- [ ] Ordenei valores numéricos com `sort -n`

---

## Parte 3 — Trabalhando com Loops

### 3.1 Atribuição de variáveis e loops `for`

```bash
# Atribuição de variáveis simples
frutas="banana abacaxi laranja"

for fruta in $frutas; do
    echo "A fruta é: $fruta"
done

# A cada iteração, "fruta" recebe um dos valores de "$frutas"
# O loop termina quando todos os valores forem processados
```

---

### 3.2 Variáveis numéricas com `declare -i`

```bash
declare -i contador=1    # Declara variável do tipo inteiro
contador=contador+1      # Incrementa o valor
echo $contador
```

```bash
declare -i i=1

while [ $i -le 3 ]; do
    echo "Iteração $i"
    i=i+1
done
```

**Operadores de comparação numérica:**

| Operador | Significado |
|----------|-------------|
| `-eq` | Igual a |
| `-ne` | Diferente de |
| `-lt` | Menor que |
| `-le` | Menor ou igual a |
| `-gt` | Maior que |
| `-ge` | Maior ou igual a |

---

### 3.3 Variáveis do tipo lista (arrays)

```bash
declare -a especies
especies=(Homo_sapiens Mus_musculus Canis_lupus)

echo ${especies[0]}   # Acessa o primeiro elemento (índice 0)
echo ${especies[1]}   # Acessa o segundo elemento (índice 1)
```

```bash
for especie in "${especies[@]}"; do
    echo "Espécie: $especie"
done

# "${especies[@]}" expande todos os elementos do array
# As aspas preservam elementos com espaços como unidades independentes
```

---

### 3.4 Trabalhando com arquivos FASTA

```bash
cd ~/2.trabalhando_com_loops/2.1_fastas
```

```bash
# Loop sobre todos os arquivos .faa no diretório
for i in *faa; do
    identificador=$(basename "$i" .faa)   # Remove a extensão .faa
    identificadores+=($identificador)     # Adiciona ao array
    grep ">" $i > ${identificador}.list  # Extrai cabeçalhos FASTA
done
```

---

### 3.5 Contando o número de sequências

```bash
grep -c ">" *   # Conta cabeçalhos FASTA em todos os arquivos
```

---

### 3.6 Verificando nomes ausentes entre arquivos

```bash
# 1. Criando uma lista de todos os IDs
grep ">" *list | cut -d ":" -f2 | sed "s/_artdb//g" | sort | uniq > lista_completa.list

# 2. Verificando IDs ausentes em um arquivo específico
grep -v -f gene6.list lista_completa.list

# 3. Formatando os dados
sed -i "s/_artdb//g" *list

# 4. Implementando em loop para todas as amostras
for i in ${identificadores[@]}; do
    grep -v -f ${i}.list lista_completa.list > ${i}.lack
done
```

**Flags do `grep` usadas aqui:**

| Flag | Função |
|------|--------|
| `-f arquivo` | Lê os padrões de busca a partir de um arquivo |
| `-v` | Inverte a lógica: seleciona linhas que NÃO correspondem |
| `-c` | Conta o número de linhas correspondentes |

> **Resultado:** O arquivo `<identificador>.lack` contém os IDs presentes em `lista_completa.list` mas **ausentes** em `<identificador>.list`.

### Checklist da Parte 3

- [ ] Criei e utilizei variáveis simples e numéricas
- [ ] Escrevi um loop `for` sobre uma lista de valores
- [ ] Escrevi um loop `while` com condição numérica
- [ ] Declarei e acessei elementos de um array
- [ ] Processei arquivos FASTA em batch com um loop
- [ ] Identifiquei IDs ausentes com `grep -v -f`

---

## Parte 4 — Testando Condições

### 4.1 Estruturas condicionais com `if`

```bash
if [ 1 -eq 1 ]; then
    echo "Isso é verdade"
fi
```

```bash
declare -i valor=10

if [[ $valor -lt 5 ]]; then
    echo "menor que 5"
elif [[ $valor -eq 10 ]]; then
    echo "igual a 10"
fi
```

**Estrutura geral do `if`:**

```bash
if [[ condição ]]; then
    # comandos se verdadeiro
elif [[ outra condição ]]; then
    # comandos se a segunda condição for verdadeira
else
    # comandos se todas as condições forem falsas
fi
```

---

### 4.2 `if` aplicado à leitura de arquivo tabular (BLAST)

```bash
cd ~/3.testando_condicoes/3.1_blast/
```

```bash
while read qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore;
    do if [[ "$pident" == "100.000" && "$length" -gt 100 ]];
    then echo "$sseqid"; fi;
    done < arquivo3.result
```

**Colunas do formato tabular BLAST (-outfmt 6):**

| Variável | Significado |
|----------|-------------|
| `qseqid` | ID da sequência de consulta (query) |
| `sseqid` | ID da sequência do banco (subject) |
| `pident` | Porcentagem de identidade |
| `length` | Comprimento do alinhamento |
| `mismatch` | Número de mismatches |
| `gapopen` | Número de aberturas de gap |
| `qstart/qend` | Início e fim na query |
| `sstart/send` | Início e fim no subject |
| `evalue` | Valor estatístico do alinhamento |
| `bitscore` | Score do alinhamento |

**Condições avaliadas:**

```
"$pident" == "100.000"   →  identidade de 100%
"$length" -gt 100        →  alinhamento com mais de 100 bases
```

Apenas quando **ambas** são verdadeiras, o `sseqid` é impresso.

---

### 4.3 Estrutura condicional `case`

```bash
opcao="A"

case $opcao in
    A)
        echo "Você escolheu A"
        ;;
    B)
        echo "Você escolheu B"
        ;;
    *)
        echo "Opção desconhecida"
        ;;
esac
```

> O `*)` funciona como "qualquer outro valor" — equivalente ao `else` no `if`.

---

### 4.4 `case` aplicado para checagem de espécies em FASTA

```bash
cd ~/3.testando_condicoes/3.2_fastas/
```

```bash
while read species; do

    case $(grep -c "$species" arquivo4.fasta) in
        0)
            echo "${species} ausente no arquivo"
            ;;
        *)
            echo "${species} presente no arquivo"
            ;;
    esac

done < generos
```

**Fluxo de execução:**

1. Lê uma espécie por vez do arquivo `generos`
2. Conta quantas vezes ela aparece em `arquivo4.fasta` com `grep -c`
3. O `case` avalia o número retornado:
   - `0` → espécie ausente
   - qualquer outro valor → espécie presente

### Checklist da Parte 4

- [ ] Escrevi um `if` simples com operador de comparação numérica
- [ ] Usei `elif` e `else` para múltiplas condições
- [ ] Processei um resultado tabular do BLAST com `while read`
- [ ] Filtrei resultados com `&&` (E lógico) no `if`
- [ ] Usei `case` para avaliar múltiplos valores possíveis
- [ ] Combinei `while read` com `case` para checar espécies em FASTA

---

## :material-flag-checkered: Você chegou ao fim do tutorial 🎉

Parabéns! Você acaba de percorrer os fundamentos do Bash aplicados à bioinformática — desde a navegação no servidor até automação com loops e condicionais.

```mermaid
flowchart LR
    A[Navegação<br/>pwd / ls / cd] --> B[Manipulação<br/>cp / mv / rm]
    B --> C[Visualização<br/>head / tail / grep]
    C --> D[Pipelines<br/>cut / sed / sort]
    D --> E[Loops<br/>for / while]
    E --> F[Condicionais<br/>if / case]
    F --> G[Automação<br/>de análises]
```

### Para se aprofundar

- **AWK** — linguagem de processamento de texto poderosa para tabelas e dados estruturados
- **Scripts Bash** — organizar pipelines em arquivos `.sh` reutilizáveis com argumentos
- **Conda/Mamba** — gerenciamento de ambientes e ferramentas bioinformáticas
- **Nextflow/Snakemake** — sistemas de workflow para pipelines reproduzíveis em larga escala
- **HPC/SLURM** — submissão de jobs em clusters de alto desempenho