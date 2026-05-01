# Filogenia molecular: inferência evolutiva com dados moleculares

!!! abstract "Resumo"
    Tutorial de inferência filogenética a partir de sequências moleculares, com aplicação em **atribuição taxonômica** — identificar a que espécie pertencem sequências "misteriosas" comparando-as com referências conhecidas.
    Cobrimos alinhamento, filtragem, concatenação, escolha de modelo evolutivo e visualização da árvore.

## :material-target: Objetivos de aprendizagem

Ao final deste tutorial, você será capaz de:

- [ ] Alinhar sequências de DNA e proteína com MAFFT e MUSCLE
- [ ] Fazer alinhamento por códon usando pal2nal
- [ ] Filtrar regiões de baixa qualidade com trimAL
- [ ] Concatenar alinhamentos de múltiplos genes
- [ ] Inferir uma árvore filogenética com IQ-TREE 3 (com e sem outgroup)
- [ ] Usar partições para genes que evoluem em taxas diferentes
- [ ] Identificar o modelo evolutivo ótimo via ModelFinderPlus
- [ ] Visualizar a árvore final no iTOL

## :material-clock-outline: Carga horária 

**4 horas**

## :material-tools: Pré-requisitos

| Item                | Detalhes                                                                         |
| ------------------- | -------------------------------------------------------------------------------- |
| Curso anterior      | [Bash / Linux para bioinfo](../bash/index.md)                                    |
| Biologia molecular  | Conceitos básicos sobre genes mitocondriais (COI, CYTB), códon e tradução        |
| Softwares           | MAFFT, MUSCLE, pal2nal, trimAL, IQ-TREE 3, catfasta2phyml                        |
| Visualização        | [Jalview](https://www.jalview.org/), [iTOL](https://itol.embl.de/)               |
| Dados               | Arquivos FASTA disponíveis em `data/exemplos/filogenia/` deste repositório       |

### Como instalar os softwares

=== "Conda (recomendado)"
```bash
    conda create -n filogenia -c bioconda mafft muscle trimal iqtree pal2nal
    conda activate filogenia
```

=== "Ubuntu/Debian (apt)"
```bash
    sudo apt-get install mafft muscle trimal iqtree
    # pal2nal e catfasta2phyml: instalação manual
```

## :material-database: Dados do exercício

Os arquivos do exercício estão na pasta [`data/exemplos/filogenia/`](https://github.com/LGBIO-UFG/PRO-BIOINFO/tree/main/data/exemplos/filogenia) deste repositório.

### Opção A — Clonar o repositório inteiro (recomendado)

```bash
git clone https://github.com/LGBIO-UFG/PRO-BIOINFO.git
cd PRO-BIOINFO/data/exemplos/filogenia/
```

### Opção B — Baixar só os arquivos do exercício

```bash
mkdir -p ex1 && cd ex1
wget https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/COI.faa
wget https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/COI.fna
wget https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/CYTB.faa
wget https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/CYTB.fna
wget https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/COI-mist.fna
wget https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/CYTB-mist.fna
```

### Conteúdo dos arquivos

| Arquivo | Conteúdo | Download |
| ------- | -------- | -------- |
| `COI.faa`         | COI referências (aa)        | [:material-download:](https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/COI.faa) |
| `COI.fna`         | COI referências (nt)        | [:material-download:](https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/COI.fna) |
| `CYTB.faa`        | CYTB referências (aa)       | [:material-download:](https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/CYTB.faa) |
| `CYTB.fna`        | CYTB referências (nt)       | [:material-download:](https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/CYTB.fna) |
| `COI-mist.fna`    | COI misteriosas (nt)        | [:material-download:](https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/COI-mist.fna) |
| `CYTB-mist.fna`   | CYTB misteriosas (nt)       | [:material-download:](https://raw.githubusercontent.com/LGBIO-UFG/PRO-BIOINFO/main/data/exemplos/filogenia/CYTB-mist.fna) |

!!! info "Cenário do exercício"
    Vamos imaginar que você recebeu duas sequências de origem desconhecida (uma de COI, outra de CYTB) e quer descobrir a que espécie pertencem, comparando-as com sequências de mamíferos conhecidos.
    

---

## :material-numeric-1-circle: Etapa 1 — Alinhamento das sequências

Antes de qualquer inferência filogenética, precisamos **alinhar as sequências** — colocar bases ou aminoácidos homólogos na mesma coluna. Esse passo é o que permite, depois, comparar posição por posição e estimar quantas mudanças evolutivas separam cada par de sequências.

!!! info "Por que alinhar antes de inferir a árvore?"
    O alinhamento garante que posições homólogas — isto é, que derivam do mesmo ancestral — fiquem alinhadas na mesma coluna. Gaps representam eventos de inserção ou deleção ao longo da evolução. Como os métodos filogenéticos comparam posição a posição, eles assumem que todas as sequências têm o mesmo comprimento e correspondência entre sítios, o que só é possível após um alinhamento adequado.

### 1.1 Conhecendo os arquivos de entrada

Os dados do exercício incluem dois genes mitocondriais (**COI** e **CYTB**) das espécies de referência, em duas formas — proteína e nucleotídeo — mais as sequências misteriosas a serem identificadas:

| Arquivo            | Conteúdo                                       | Quando usar                                  |
| ------------------ | ---------------------------------------------- | -------------------------------------------- |
| `COI.faa`          | COI das referências, em **aminoácidos**        | Alinhamento de proteína                      |
| `COI.fna`          | COI das referências, em **nucleotídeos**       | Alinhamento de nucleotídeo                   |
| `CYTB.faa`         | CYTB das referências, em aminoácidos           | Alinhamento de proteína                      |
| `CYTB.fna`         | CYTB das referências, em nucleotídeos          | Alinhamento de nucleotídeo                   |
| `COI-mist.fna`     | COI das sequências **misteriosas** (DNA)       | Comparação contra o conjunto de referência   |
| `CYTB-mist.fna`    | CYTB das sequências misteriosas (DNA)          | Comparação contra o conjunto de referência   |

!!! tip "Proteína ou nucleotídeo: qual escolher?"
    Alinhamentos de proteína são geralmente mais confiáveis para espécies distantes, pois os aminoácidos evoluem mais lentamente e preservam melhor a homologia. Já os nucleotídeos carregam mais variação, sendo úteis para distinguir espécies próximas. O alinhamento por códon combina essas vantagens: usa a robustez da proteína para alinhar e mantém a informação completa do DNA.

### 1.2 Alinhamento de proteína

Vamos começar pelo gene COI em sua versão de proteína (`COI.faa`). Existem várias ferramentas — duas das mais usadas são MAFFT e MUSCLE:

=== "MAFFT (recomendado)"
    ````bash
    mafft COI.faa > COI.ali.faa
    ````

    O MAFFT é rápido e tem várias estratégias automáticas. Sem flags, ele escolhe sozinho a estratégia mais adequada ao tamanho do conjunto.

=== "MUSCLE"
    ```bash
    muscle -in COI.faa -out COI.muscle.ali.faa
    ```

    O MUSCLE é outro alinhador clássico. O MUSCLE pode produzir alinhamentos ligeiramente mais precisos em conjuntos pequenos e bem conservados. Já o MAFFT é mais rápido e escala melhor para conjuntos grandes ou mais divergentes, sendo geralmente a escolha padrão em pipelines modernos.

!!! info "Convenção de nomes"
    Usamos o sufixo `.ali` para deixar claro que é a versão **alinhada** do arquivo. Convenções como essa salvam horas no futuro: você sabe o que é cada arquivo só pelo nome.

### 1.3 Alinhamento de nucleotídeo

Para alinhar diretamente as sequências de nucleotídeo do COI:

```bash
mafft COI.fna > COI.ali.fna
```

O alinhamento direto de nucleotídeos é adequado para sequências próximas, onde a estrutura de códons ainda está preservada. Em conjuntos mais divergentes, porém, substituições silenciosas e indels podem desalinhar códons, levando a inferências incorretas — nesses casos, o alinhamento por proteína ou por códon é mais robusto.

### 1.4 Alinhamento por códon (pal2nal)

Quando temos tanto a proteína quanto o nucleotídeo da **mesma região**, podemos fazer um alinhamento **por códon**: alinhar pelas proteínas (mais conservadas) e depois projetar de volta nos nucleotídeos. Isso preserva a estrutura de leitura (ORF) e dá um alinhamento de DNA muito mais limpo.

A ferramenta pal2nal faz exatamente isso:

```bash
pal2nal COI.ali.faa COI.fna -output fasta -codontable 2 > COI.ali.codon.fa
```

!!! warning "Atenção: tabela de códon mitocondrial"
    O parâmetro `-codontable 2` é **crítico** — usa a **tabela genética mitocondrial de vertebrados**, que difere da tabela padrão (a número 1) em alguns códons importantes:

    Por exemplo, o códon TGA, que na tabela padrão codifica um sinal de parada, na mitocôndria de vertebrados codifica o aminoácido triptofano (W). Já os códons AGA e AGG, que normalmente codificam arginina, funcionam como códons de parada na mitocôndria. Essas diferenças tornam essencial o uso da tabela correta.

    Se você esquecer essa flag, o pal2nal pode interpretar códons como _stop_ e abortar a tradução no meio.

??? info "Outras tabelas de códon comuns"
    Lista completa: [Genetic Codes do NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

    Resumo das mais usadas em bioinformática:

    | Número | Tabela                                       | Uso típico                       |
    | ------ | -------------------------------------------- | -------------------------------- |
    | 1      | Standard Code                                | Genoma nuclear da maioria        |
    | 2      | Vertebrate Mitochondrial                     | Mitocôndria de vertebrados       |
    | 4      | Mold/Protozoan/Coelenterate Mitochondrial    | Mitocôndria de fungos, parasitas |
    | 5      | Invertebrate Mitochondrial                   | Mitocôndria de invertebrados     |
    | 11     | Bacterial, Archaeal and Plant Plastid        | Bactérias, arqueias, plastídios  |

### 1.5 Visualizando o alinhamento

Antes de seguir, vale **olhar** o alinhamento — alinhamento ruim no início = árvore errada no fim. Uma das opções mais práticas é o [Jalview](https://www.jalview.org/), que abre arquivos FASTA alinhados e mostra colunas, conservação, qualidade.

!!! tip "O que olhar no Jalview"
    - **Colunas bem conservadas** (cores uniformes) indicam regiões homólogas confiáveis.
    - **Regiões com muitos gaps espalhados** podem ser indicadores de erro de alinhamento — vale considerar trimming (próximo passo) ou re-alinhar com outros parâmetros.
    - **Sequências com muito gap** isoladas das outras podem ser parálogas, contaminação ou de baixa qualidade — considere remover.

### :material-pencil-outline: Exercício 1.1

Alinhe **também** o gene CYTB nas três modalidades (proteína com MAFFT, nucleotídeo com MAFFT, e por códon com pal2nal). Os arquivos de saída devem se chamar:

- `CYTB.ali.faa`
- `CYTB.ali.fna`
- `CYTB.ali.codon.fa`

??? success "Comandos esperados"
    ```bash
    # Proteína com MAFFT
    mafft CYTB.faa > CYTB.ali.faa

    # Nucleotídeo com MAFFT
    mafft CYTB.fna > CYTB.ali.fna

    # Códon com pal2nal (atenção à tabela 2)
    pal2nal CYTB.ali.faa CYTB.fna -output fasta -codontable 2 > CYTB.ali.codon.fa
    ```

### Checklist da Etapa 1

- [ ] Tenho `COI.ali.faa`, `COI.ali.fna`, `COI.ali.codon.fa`
- [ ] Tenho `CYTB.ali.faa`, `CYTB.ali.fna`, `CYTB.ali.codon.fa`
- [ ] Inspecionei pelo menos um dos alinhamentos no Jalview

---

## :material-numeric-2-circle: Etapa 2 — Alinhar misteriosas com as referências

Para que a inferência filogenética posicione as sequências misteriosas dentro da árvore das referências, **ambas precisam estar no mesmo alinhamento múltiplo**. Alinhar separadamente quebraria a homologia posicional.

!!! info "Por que alinhar tudo junto?"
    Cada coluna do alinhamento representa uma posição homóloga compartilhada entre todos os táxons analisados. Se as sequências misteriosas e as referências fossem alinhadas separadamente, as colunas de cada alinhamento não corresponderiam necessariamente aos mesmos sítios evolutivos — o que comprometeria a comparabilidade e introduziria erro sistemático na inferência filogenética.

### 2.1 Combinar referências + misteriosas

Os arquivos do exercício já vêm com referências e misteriosas no mesmo FASTA: `COI-mist.fna` e `CYTB-mist.fna`. Por isso podemos pular direto para o alinhamento.

!!! tip "Se as sequências estiverem separadas em outros projetos"
    Em casos reais, suas misteriosas podem estar em um arquivo `mist.fna` à parte das referências. Para juntá-las antes de alinhar, use `cat`:

```bash
    cat COI.fna mist.fna > COI-mist.fna
```

    O `cat` concatena arquivos em sequência. Funciona com FASTA porque cada sequência é um bloco independente começando com `>` — basta colocar todos os blocos no mesmo arquivo e o MAFFT lê normalmente.

### 2.2 Alinhar com MAFFT

```bash
# Gene COI (referências + misteriosas)
mafft COI-mist.fna > COI-mist.ali.fna

# Gene CYTB (referências + misteriosas)
mafft CYTB-mist.fna > CYTB-mist.ali.fna
```

Aqui usamos o MAFFT diretamente sobre os nucleotídeos. Como o conjunto envolve mamíferos no mesmo gene mitocondrial (espécies relativamente próximas), o alinhamento de nucleotídeo direto já é confiável — não precisamos da etapa intermediária com `pal2nal` aqui.

### 2.3 Confirmação rápida no Jalview

Antes de seguir para a filtragem, abre os dois alinhamentos no [Jalview](https://www.jalview.org/) e confere:

- As **misteriosas** se alinharam dentro das mesmas colunas conservadas das referências?
- Há regiões com muito gap no início ou no fim, indicando que ali precisa ser cortado na próxima etapa?
- Alguma misteriosa está **muito divergente** das referências (gaps em massa)? Se sim, vale anotar — pode ser sequência distante, contaminação ou parálogo.

### :material-pencil-outline: Exercício 2.1

Quantas sequências há em cada arquivo combinado (referências + misteriosas)? Use `grep` para contar.

??? success "Resposta"
    ```bash
    grep -c "^>" COI-mist.fna
    grep -c "^>" CYTB-mist.fna
    ```

    Cada `^>` corresponde ao cabeçalho de uma sequência FASTA.

### Checklist da Etapa 2

- [ ] Tenho `COI-mist.ali.fna`
- [ ] Tenho `CYTB-mist.ali.fna`
- [ ] Inspecionei pelo menos um dos alinhamentos no Jalview
- [ ] Confirmei que as misteriosas se alinham bem entre as referências

---

## :material-numeric-3-circle: Etapa 3 — Filtragem do alinhamento (trimAL)

Nem toda coluna do alinhamento contribui igualmente para a inferência da árvore. Regiões com muitos gaps, baixa conservação ou ruído aleatório podem **atrapalhar mais do que ajudar** — o sinal filogenético fica diluído. A filtragem (também chamada de _trimming_) remove essas colunas problemáticas antes de seguir para a inferência.

!!! info "Por que filtrar?"
    Colunas com muitos gaps ou alta variabilidade não informativa tendem a introduzir ruído, dificultando a recuperação do sinal evolutivo real. Ao remover essas posições, aumentamos a proporção de sítios informativos, o que geralmente resulta em árvores mais estáveis e com melhor suporte estatístico.

### 3.1 trimAL com método automático

O **trimAL** é uma das ferramentas mais usadas para esse fim. A flag `-automated1` aplica uma heurística que escolhe automaticamente entre `gappyout` e `strict` com base nas características do alinhamento — funciona bem como ponto de partida em quase todos os casos.

```bash
# Gene COI
trimal -in COI-mist.ali.fna -automated1 -out COI-mist.trim.ali.fna

# Gene CYTB
trimal -in CYTB-mist.ali.fna -automated1 -out CYTB-mist.trim.ali.fna
```

??? info "Outros métodos do trimAL (clique para expandir)"
    O `-automated1` é conveniente, mas vale conhecer alternativas para casos específicos:

    | Flag             | O que faz                                                                |
    | ---------------- | ------------------------------------------------------------------------ |
    | `-gappyout`      | Remove colunas com muitos gaps; bom para alinhamentos com muitas inserções |
    | `-strict`        | Critério mais agressivo; bom quando o alinhamento já é muito conservado   |
    | `-strictplus`    | Variação ainda mais conservadora do strict                               |
    | `-gt 0.8`        | Mantém apenas colunas com no máximo 20% de gaps (cutoff manual)          |
    | `-cons 60`       | Garante que pelo menos 60% das colunas sejam mantidas                    |

    Documentação completa: [trimAL manual](https://vicfero.github.io/trimal/).

### 3.2 Inspeção visual após o trimming

Abre os arquivos `*.trim.ali.fna` no Jalview e compara lado a lado com os originais (`*.ali.fna`). Você deve ver:

- **Menos colunas** no total (regiões problemáticas foram removidas)
- **Maior consistência** nas colunas mantidas (cores mais uniformes nas posições de cada base)
- As **misteriosas continuam alinhadas** dentro das referências, agora em colunas mais limpas

!!! warning "Atenção: trimAL pode ser agressivo demais"
    Em alinhamentos curtos ou com poucas sequências, o `-automated1` pode remover mais colunas do que o desejado. Sempre verifique se o alinhamento final ainda tem comprimento suficiente (idealmente algumas centenas de posições) e se as sequências misteriosas continuam bem representadas — caso contrário, prefira um critério mais permissivo.

### :material-pencil-outline: Exercício 3.1

Quantas colunas tinha o alinhamento de COI **antes** e **depois** da filtragem?

??? success "Comando esperado"
    ```bash
    # Imprime o tamanho da primeira sequência = número de colunas do alinhamento
    awk '/^>/{if(seq){print length(seq); exit} seq=""; next} {seq=seq $0}' COI-mist.ali.fna
    awk '/^>/{if(seq){print length(seq); exit} seq=""; next} {seq=seq $0}' COI-mist.trim.ali.fna
    ```

    Como todas as sequências de um alinhamento têm o mesmo comprimento (incluindo gaps), basta medir uma. O segundo número deve ser **menor** — quanto, depende do quão "sujo" estava o alinhamento original.

### Checklist da Etapa 3

- [ ] Tenho `COI-mist.trim.ali.fna`
- [ ] Tenho `CYTB-mist.trim.ali.fna`
- [ ] Comparei pelo menos um deles antes/depois no Jalview
- [ ] O alinhamento filtrado tem colunas suficientes para inferir a árvore

---

## :material-numeric-4-circle: Etapa 4 — Concatenar os alinhamentos dos dois genes

Quando temos múltiplos genes alinhados separadamente, juntá-los em um **único super-alinhamento** (chamado de **supermatriz**) costuma resultar em uma árvore mais robusta. A ideia é simples: mais dados por táxon → mais sinal filogenético → resolução melhor de relações que cada gene sozinho talvez não conseguisse distinguir.

!!! info "Por que concatenar?"
    Genes diferentes podem evoluir em taxas distintas, mas compartilham, em geral, a mesma história evolutiva das espécies. Ao concatená-los, aumentamos a quantidade de informação por táxon, o que tende a melhorar a resolução da topologia e o suporte estatístico da árvore — especialmente quando cada gene isolado tem sinal limitado.

### 4.1 Concatenando com catfasta2phyml

A ferramenta `catfasta2phyml.pl` faz a concatenação respeitando os táxons: combina sequências do **mesmo táxon** entre os arquivos e preenche com gaps quando algum táxon falta em algum gene.

```bash
catfasta2phyml.pl -c -f *mist.trim.ali.fna > concatenado.fasta
```

O que cada pedaço significa:

- `-c` — modo de **concatenação** (sem isso, o script só converte formato)
- `-f` — saída em formato **FASTA** (sem essa flag, sai em Phylip)
- `*mist.trim.ali.fna` — wildcard que pega `COI-mist.trim.ali.fna` **e** `CYTB-mist.trim.ali.fna`
- `> concatenado.fasta` — redireciona a saída para um arquivo

!!! warning "⚠️ ANOTE as posições dos genes que aparecerem no terminal"
    Depois de rodar o comando, o `catfasta2phyml.pl` imprime na tela algo parecido com:

    ```text
    COI-mist.trim.ali.fna     = 1-581
    CYTB-mist.trim.ali.fna    = 582-1281
    ```

    Esses números — **posições inicial e final** de cada gene no super-alinhamento — vão ser usados na **Etapa 6** (filogenia com partições). **Anote agora** num arquivo `partitions-info.txt` ou em qualquer lugar que você não vá perder. Sem essas posições, não dá pra montar o arquivo de partição depois.

    !!! tip "Truque: salvar a saída automaticamente"
        Pra não depender de copiar manualmente, redirecione também o STDERR:

    ```bash
        catfasta2phyml.pl -c -f *mist.trim.ali.fna > concatenado.fasta 2> partitions-info.txt
    ```

        O `2>` captura as mensagens de erro/info (que é onde o script imprime as posições) num arquivo separado.

### 4.2 Verificação rápida

```bash
# Quantas sequências há no concatenado?
grep -c "^>" concatenado.fasta

# Quantas colunas tem o super-alinhamento?
awk '/^>/{if(seq){print length(seq); exit} seq=""; next} {seq=seq $0}' concatenado.fasta
```

- O número de **sequências** deve bater com o de qualquer um dos genes individuais (cada táxon vira uma linha do super-alinhamento).
- O número de **colunas** deve ser a **soma** das colunas dos dois alinhamentos filtrados da Etapa 3.

!!! info "E se algum táxon faltar em um dos genes?"
    O `catfasta2phyml.pl` preenche automaticamente com gaps a região correspondente ao gene ausente para aquele táxon. Isso é aceitável até certo ponto, mas uma alta proporção de dados ausentes (por exemplo, >50% do alinhamento) pode reduzir a confiabilidade da inferência — nesses casos, vale considerar remover o táxon ou analisar os genes separadamente.

### :material-pencil-outline: Exercício 4.1

Confirme que a concatenação fez sentido: número de táxons igual ao dos arquivos individuais, e número de colunas igual à soma.

??? success "Comandos esperados"
    ```bash
    # Sequências em cada arquivo
    echo "COI:          $(grep -c '^>' COI-mist.trim.ali.fna)"
    echo "CYTB:         $(grep -c '^>' CYTB-mist.trim.ali.fna)"
    echo "Concatenado:  $(grep -c '^>' concatenado.fasta)"

    # Colunas em cada arquivo
    awk '/^>/{if(seq){print "COI:         " length(seq); exit} seq=""; next} {seq=seq $0}' COI-mist.trim.ali.fna
    awk '/^>/{if(seq){print "CYTB:        " length(seq); exit} seq=""; next} {seq=seq $0}' CYTB-mist.trim.ali.fna
    awk '/^>/{if(seq){print "Concat:      " length(seq); exit} seq=""; next} {seq=seq $0}' concatenado.fasta
    ```

    Esperado: as três contagens de sequência são iguais, e as colunas do concatenado = COI + CYTB.

### Checklist da Etapa 4

- [ ] Tenho `concatenado.fasta`
- [ ] **Anotei** as posições (inicial-final) de cada gene em `partitions-info.txt` ou similar
- [ ] Conferi que o número de táxons e o comprimento batem

---

## :material-numeric-5-circle: Etapa 5 — Inferência da árvore filogenética (IQ-TREE 3)

Com o super-alinhamento pronto, finalmente podemos **inferir a árvore**. Vamos usar o **IQ-TREE 3**, que implementa Máxima Verossimilhança (ML) e é hoje uma das principais opções em inferência filogenética molecular — combina acurácia, velocidade e seleção automática de modelo.

!!! info "O que o IQ-TREE está fazendo?"
    O IQ-TREE busca, entre as possíveis árvores, aquela que melhor explica os dados observados sob um modelo evolutivo — isto é, a árvore com maior verossimilhança. Durante esse processo, ele também estima os parâmetros do modelo e avalia o suporte estatístico dos clados por meio de reamostragem (bootstrap). Como o número de árvores possíveis cresce rapidamente com o número de táxons, essa busca é feita por métodos heurísticos, explorando o espaço de soluções de forma eficiente.

### 5.1 Inferência sem outgroup

A versão mais simples roda o IQ-TREE sem especificar nenhum táxon como referência externa:

```bash
iqtree3 -s concatenado.fasta -m MFP -pre arvore-ex1-sem-out -B 1000
```

| Flag                       | Significado                                                              |
| -------------------------- | ------------------------------------------------------------------------ |
| `-s concatenado.fasta`     | Arquivo de entrada (super-alinhamento)                                   |
| `-m MFP`                   | **ModelFinder Plus** — escolhe automaticamente o melhor modelo evolutivo |
| `-pre arvore-ex1-sem-out`  | Prefixo dos arquivos de saída                                            |
| `-B 1000`                  | **Ultrafast Bootstrap** com 1000 réplicas para suporte de nós            |

!!! tip "Arquivos gerados pelo IQ-TREE"
    O comando produz vários arquivos com o prefixo escolhido. Os principais:

    | Arquivo            | Conteúdo                                              |
    | ------------------ | ----------------------------------------------------- |
    | `*.treefile`       | A árvore final em formato Newick (texto)              |
    | `*.iqtree`         | Relatório completo (modelo, ML score, suportes)       |
    | `*.log`            | Log da execução                                       |
    | `*.contree`        | Árvore de consenso construída a partir do bootstrap   |
    | `*.bionj`          | Árvore inicial (BioNJ) usada como semente             |

### 5.2 Inferência com outgroup

Quando temos um táxon que sabemos ser **mais distante** do grupo de interesse, podemos usá-lo como **outgroup** para enraizar a árvore. Neste exercício, a capivara (*Hydrochoerus hydrochaeris*) faz esse papel:

```bash
iqtree3 -s concatenado.fasta -m MFP -o Hydrochoerus_hydrochaeris -pre arvore-ex1-com-out -B 1000
```

A única flag adicional é `-o`:

| Flag                              | Significado                                                  |
| --------------------------------- | ------------------------------------------------------------ |
| `-o Hydrochoerus_hydrochaeris`    | Táxon de outgroup (use o **nome exato** do cabeçalho FASTA, sem espaços) |

!!! info "Para que serve o outgroup?"
    O outgroup permite enraizar a árvore, definindo onde está a raiz e, portanto, a direção das divergências evolutivas. Sem um outgroup, a árvore é não-enraizada: ela mostra as relações entre os táxons, mas não indica qual linhagem se separou primeiro.

!!! warning "Como escolher um bom outgroup"
    Um bom outgroup deve ser filogeneticamente próximo o suficiente para permitir um alinhamento confiável, mas distante o bastante para não pertencer ao grupo de interesse. Se for muito distante, pode introduzir erros por saturação ou alinhamento ruim; se for muito próximo, pode acabar sendo agrupado dentro do ingroup. Idealmente, deve haver evidência prévia de que ele é externo ao grupo analisado.
    
### 5.3 Interpretando os valores de suporte (bootstrap)

Cada nó da árvore final vem rotulado com um número de **0 a 100** — a porcentagem de réplicas de bootstrap que recuperaram aquele clado.

| Valor de suporte    | Interpretação                                                  |
| ------------------- | -------------------------------------------------------------- |
| **≥ 95**            | Clado bem suportado — alta confiança                           |
| **80 – 95**         | Suporte moderado, mas razoável                                 |
| **50 – 80**         | Suporte fraco — relação dúbia, considere ambígua               |
| **< 50**            | Efetivamente sem suporte                                       |

!!! tip "Atenção: ultrafast bootstrap tende a inflar valores"
    O `-B` (ultrafast bootstrap) é mais rápido que o bootstrap clássico (`-b`), mas **gera valores mais altos** para o mesmo dado. A regra de bolso é considerar **≥ 95** como bom suporte com `-B`, em vez dos ≥ 70 que seriam o cutoff clássico.

### :material-pencil-outline: Exercício 5.1

Qual foi o modelo evolutivo escolhido automaticamente pelo ModelFinderPlus?

??? success "Comando esperado"
    ```bash
    grep -i "best-fit" arvore-ex1-com-out.log
    ```

    A saída mostra algo do tipo `Best-fit model: GTR+F+I+G4 chosen according to BIC`. Esse modelo é o que o IQ-TREE usou em toda a inferência. Vamos voltar nele na Etapa 7.

### Checklist da Etapa 5

- [ ] Tenho `arvore-ex1-sem-out.treefile` (sem outgroup)
- [ ] Tenho `arvore-ex1-com-out.treefile` (com outgroup)
- [ ] Identifiquei o modelo evolutivo escolhido pelo ModelFinder
- [ ] Conferi que rodou sem erro (arquivo `*.iqtree` foi gerado)

---

## :material-numeric-6-circle: Etapa 6 — Filogenia com partições por gene

Quando concatenamos dois genes em um único super-alinhamento (Etapa 4), aplicar **um único modelo evolutivo** sobre toda a matriz pressupõe que ambos evoluem da mesma forma — o que raramente é verdade. Genes mitocondriais como COI e CYTB podem ter taxas, vieses de composição e padrões de substituição diferentes, por causa de pressões seletivas e papéis funcionais distintos.

A solução é dividir o alinhamento em **partições** (uma por gene) e deixar o IQ-TREE estimar **um modelo por partição**.

!!! info "Por que usar partições?"
    Cada gene pode ter sua própria taxa evolutiva, composição de bases e padrões de substituição, refletindo diferentes pressões seletivas e funções biológicas. Aplicar um único modelo a todo o alinhamento ignora essas diferenças e pode distorcer a inferência. Ao usar partições, cada gene é modelado de forma independente, e a árvore resultante integra essas informações de maneira mais realista e biologicamente consistente.

### 6.1 Criar o arquivo de partições

O IQ-TREE entende partições no formato **Nexus**. Você precisa montar um arquivo `partitions.nex` no formato abaixo, **substituindo os números** pelas posições que você anotou na Etapa 4:

```text
#nexus
begin sets;
  charset COI  = 1-581;
  charset CYTB = 582-1281;
end;
```

O que esse arquivo está dizendo:

- `#nexus` — declaração obrigatória do formato
- `begin sets; ... end;` — bloco de definição de partições
- `charset NOME = INICIO-FIM` — define uma partição com nome `NOME` ocupando as colunas de `INICIO` até `FIM`

!!! warning "Use SUAS posições, não as do exemplo"
    Os números `1-581` e `582-1281` são ilustrativos. **Os seus números** vieram do `catfasta2phyml.pl` na Etapa 4 e podem ser diferentes — depende de quanto cada gene perdeu na filtragem do trimAL. Confere no seu arquivo:

    ```bash
    cat partitions-info.txt
    ```

    E coloca exatamente esses valores no `partitions.nex`.

!!! tip "Atalho: criar o arquivo direto pelo terminal"
    Em vez de abrir editor, dá pra gerar o arquivo com um here-document:

    ```bash
    cat > partitions.nex << 'EOF'
    #nexus
    begin sets;
      charset COI  = 1-581;
      charset CYTB = 582-1281;
    end;
    EOF
    ```

    O `<< 'EOF' ... EOF` é um **here-document**: tudo que estiver entre os dois `EOF` vai parar dentro do arquivo. Ajusta os números antes de rodar e confere depois com `cat partitions.nex`.

### 6.2 Rodar o IQ-TREE com partições

A única flag nova é `-p partitions.nex`:

```bash
iqtree3 -s concatenado.fasta \
        -m MFP \
        -o Hydrochoerus_hydrochaeris \
        -p partitions.nex \
        -pre arvore-ex1-com-out-e-part \
        -B 1000
```

> O `\` no fim de cada linha é só pra quebrar o comando em várias linhas e ficar legível. Em uma linha só também funciona.

A diferença em relação à Etapa 5: o IQ-TREE agora roda o ModelFinder **separadamente para cada partição**, escolhendo o melhor modelo pra COI e o melhor pra CYTB independentemente. A árvore final é estimada usando os dois modelos juntos.

### 6.3 Comparação entre as três árvores

Você agora tem **três árvores** das três rodadas do IQ-TREE:

| Arquivo                                  | Variante                                  |
| ---------------------------------------- | ----------------------------------------- |
| `arvore-ex1-sem-out.treefile`            | Sem outgroup, sem partições               |
| `arvore-ex1-com-out.treefile`            | Com outgroup, sem partições               |
| `arvore-ex1-com-out-e-part.treefile`     | Com outgroup **e** partições              |

!!! info "Qual usar para a atribuição taxonômica?"
    A árvore com outgroup e partições é a mais adequada, pois combina enraizamento correto com modelagem independente dos genes. Isso tende a produzir uma topologia mais confiável e melhor suportada, sendo a melhor base para identificar as sequências misteriosas.

### :material-pencil-outline: Exercício 6.1

O ModelFinder escolheu modelos diferentes para COI e CYTB? Compare também com o modelo único da rodada anterior.

??? success "Comandos esperados"
    ```bash
    # Modelo único (rodada da Etapa 5, com outgroup)
    grep -i "best-fit model" arvore-ex1-com-out.log

    # Modelos por partição (rodada com partições)
    grep -i "best.*model" arvore-ex1-com-out-e-part.log | head
    ```

    Na rodada com partições devem aparecer dois modelos — um para COI, outro para CYTB. Frequentemente são diferentes entre si, e podem também diferir do modelo único da rodada sem partições.

### Checklist da Etapa 6

- [ ] Criei `partitions.nex` com as **minhas** posições (não as do exemplo)
- [ ] Tenho `arvore-ex1-com-out-e-part.treefile`
- [ ] Conferi os modelos escolhidos para COI e CYTB

---

## :material-numeric-7-circle: Etapa 7 — Identificar o modelo evolutivo escolhido (ModelFinderPlus)

Em todas as rodadas usamos a flag `-m MFP` (ModelFinder Plus), que **testa dezenas de modelos** e escolhe automaticamente o melhor com base em critérios como BIC ou AIC. Vale agora olhar **qual modelo foi de fato escolhido** — tanto por curiosidade quanto pra documentar a análise.

!!! info "Por que documentar o modelo escolhido?"
    A reprodutibilidade de uma análise filogenética depende do conhecimento exato do modelo evolutivo utilizado. Reportar o modelo selecionado (e o critério de escolha, como BIC ou AIC) é uma prática essencial, tanto para permitir a reprodução dos resultados quanto para garantir transparência em análises científicas.

### 7.1 Modelo único (rodadas sem partições)

Quando rodamos **sem** `-p`, o IQ-TREE escolhe um único modelo para todo o alinhamento. Ele aparece no `*.log`:

```bash
# Para qualquer das rodadas sem partição:
grep -i "best-fit" *.log
```

Saída esperada:

```text
arvore-ex1-sem-out.log:Best-fit model: GTR+F+I+G4 chosen according to BIC
arvore-ex1-com-out.log:Best-fit model: GTR+F+I+G4 chosen according to BIC
```

### 7.2 Modelos por partição

Quando rodamos com `-p partitions.nex`, o IQ-TREE escolhe **um modelo por partição** (ou pode até sugerir agrupar partições, dependendo do caso). O esquema final é gravado num arquivo nexus separado:

```bash
more arvore-ex1-com-out-e-part.best_scheme.nex
```

!!! tip "Esse arquivo é reutilizável"
    O `.best_scheme.nex` é um nexus completo, pronto pra ser passado novamente ao IQ-TREE via `-p` numa próxima análise — pulando assim toda a etapa de seleção de modelo. Útil quando você quer rodar variações da análise (mais bootstraps, outras flags) sem refazer o ModelFinder do zero.

### 7.3 Decifrando o nome do modelo

Os nomes vêm com várias partes coladas com `+`. Os componentes mais comuns para nucleotídeos:

| Sufixo  | Significado                                                                          |
| ------- | ------------------------------------------------------------------------------------ |
| `GTR`   | **General Time Reversible** — modelo de substituição mais flexível para nucleotídeos |
| `HKY`   | Modelo mais simples (transições ≠ transversões); usado quando GTR é "exagero"        |
| `TIM2`, `TPM2`, `TVM` etc. | Variações intermediárias entre HKY e GTR                                |
| `+F`    | Frequências de bases **estimadas dos dados** (em vez de assumir 0.25 cada)           |
| `+I`    | **Proporção de sítios invariáveis** — parte do alinhamento que não muda              |
| `+G4`   | Heterogeneidade de taxa entre sítios em **4 categorias gama**                        |
| `+R3`   | Variação de taxa em 3 categorias livres (alternativa mais flexível a `+G`)           |

!!! info "Exemplo de leitura"
    `GTR+F+I+G4` quer dizer: modelo GTR, com frequências de bases estimadas dos dados (`+F`), proporção de sítios invariáveis (`+I`) e taxa heterogênea entre sítios em 4 categorias gama (`+G4`). É um dos modelos mais ricos para dados nucleotídicos e frequentemente o escolhido quando há sinal estatístico suficiente para justificar a complexidade.

### :material-pencil-outline: Exercício 7.1

Os modelos escolhidos para COI e CYTB foram iguais ou diferentes? E o que isso revela sobre os dois genes?

??? success "Comando esperado e interpretação"
    ```bash
    cat arvore-ex1-com-out-e-part.best_scheme.nex
    ```

    Você verá algo como:

    ```text
    #nexus
    begin sets;
      charset COI = 1-581;
      charset CYTB = 582-1281;
      charpartition mymodels =
        TIM2+F+G4: COI,
        HKY+F+I: CYTB;
    end;
    ```

    Modelos diferentes para COI e CYTB indicam que esses genes apresentam padrões evolutivos distintos no conjunto analisado, como diferenças em taxa de substituição ou composição de bases. Isso reforça a importância do uso de partições, pois cada gene pode ser descrito de forma mais adequada por seu próprio modelo evolutivo.

### Checklist da Etapa 7

- [ ] Identifiquei o modelo único da rodada sem partições
- [ ] Conferi o `.best_scheme.nex` da rodada com partições
- [ ] Entendi o que cada parte do nome do modelo significa

---

## :material-numeric-8-circle: Etapa 8 — Visualizar a árvore e identificar as misteriosas

Chegou o momento da verdade: visualizar a árvore filogenética e descobrir **a que espécies pertencem** as sequências misteriosas. Vamos usar o **iTOL** (Interactive Tree of Life), uma ferramenta web gratuita que renderiza árvores em formato Newick e permite muitas customizações sem instalar nada.

!!! info "Por que iTOL?"
    O iTOL é gratuito, roda diretamente no navegador e permite explorar árvores filogenéticas de forma interativa, sem necessidade de instalação. Ele suporta os principais formatos usados em filogenia e facilita tanto a visualização quanto a geração de figuras prontas para publicação.

### 8.1 Dar uma olhada na árvore como texto (Newick)

Antes de subir, vale espiar o arquivo:

```bash
cat arvore-ex1-com-out-e-part.treefile
```

Você vai ver algo do tipo:

```text
((Hydrochoerus_hydrochaeris:0.123,(Mus_musculus:0.045,Rattus_norvegicus:0.038)98:0.087)100:0.012, ... );
```

!!! info "Formato Newick em três regras"
    Cada par de parênteses representa um agrupamento (clado) de táxons relacionados. Os valores após `:` indicam o comprimento dos ramos, ou seja, a quantidade de mudança evolutiva. Já os números associados aos nós internos representam o suporte estatístico (bootstrap) daquele agrupamento.

### 8.2 Subir a árvore no iTOL

1. Acesse [https://itol.embl.de/](https://itol.embl.de/)
2. Clica em **"Upload"** no topo da página
3. Faz upload do `arvore-ex1-com-out-e-part.treefile` (ou cola o conteúdo direto)
4. Dá um nome à árvore e confirma

Em segundos a árvore aparece renderizada de forma interativa.

### 8.3 Customizações úteis no iTOL

| Ação                                  | Onde                                                       |
| ------------------------------------- | ---------------------------------------------------------- |
| Trocar o estilo (retangular, circular)| Painel direito, "Mode"                                     |
| Mostrar valores de bootstrap          | "Advanced" → "Bootstrap / display values"                  |
| Reposicionar o outgroup como raiz     | Clica no nó da capivara → "Tree structure" → "Reroot the tree here" |
| Colorir clados/ramos                  | Clica no nó → "Branch color"                              |
| Destacar as **misteriosas**           | "Datasets" → "Color labels" ou "Color ranges"              |
| Exportar (SVG/PNG/PDF)                | "Export" no painel direito                                 |

!!! tip "Destacar as misteriosas pra leitura rápida"
    Use a opção "Color labels" no iTOL para destacar as sequências misteriosas com uma cor contrastante (por exemplo, vermelho), mantendo as referências em cores neutras. Isso facilita identificar rapidamente em qual clado cada misteriosa se posiciona, especialmente em árvores maiores.

### 8.4 Interpretando: atribuição taxonômica

A pergunta-chave: **em qual clado de referência cada misteriosa se posicionou?**

Para cada misteriosa, observe:

1. **A qual clado de referência ela é irmã** (ou se está dentro do clado)
2. **O suporte de bootstrap do nó** que une a misteriosa às referências

A atribuição é robusta quando a misteriosa cai em um clado de referência conhecido com **bootstrap ≥ 95**.

!!! warning "Cuidado com clados mal suportados"
    Se a misteriosa cair em um nó com bootstrap baixo (< 80), a atribuição é incerta. **Reporte como "inconclusiva"** em vez de chutar — em projetos reais, essa honestidade evita publicar identificação errada.

!!! tip "Quando o COI e o CYTB discordam"
    Como você fez **uma única árvore concatenada**, em geral COI e CYTB já estão "votando juntos" no posicionamento da misteriosa. Mas se você quiser uma checagem extra, dá pra rodar árvores **separadas** por gene (sem concatenar) e ver se as duas concordam. Discordância forte costuma indicar paralogo, contaminação ou hibridização.

### :material-pencil-outline: Exercício 8.1

Identifique cada uma das sequências misteriosas no seu dataset. Em qual clado caíram, com qual suporte, e qual é a atribuição taxonômica?

??? success "Resposta"
    "**Misteriosa 1** caiu como irmã de *XXXX XXXX* com suporte XX% — atribuição: gênero *XXX*, possivelmente *XXX XXXX* ou espécie próxima.
    **Misteriosa 2** caiu dentro do clado de *XXXXX*, próxima a *X. XXXX*, com suporte 97% — atribuição: gênero *XXXXX*."]

### Checklist da Etapa 8

- [ ] Subi a árvore no iTOL e ela renderizou sem erro
- [ ] Reorientei a árvore com a capivara como outgroup (raiz)
- [ ] Habilitei a exibição dos valores de bootstrap
- [ ] Pintei as misteriosas para leitura rápida
- [ ] Identifiquei o clado de cada misteriosa
- [ ] Avaliei o suporte antes de fazer a atribuição
- [ ] Exportei a figura final em formato apropriado

---

## :material-flag-checkered: Você chegou ao fim do tutorial 🎉

Parabéns! Você acaba de percorrer uma **pipeline completa de inferência filogenética** para atribuição taxonômica:

```mermaid
flowchart LR
    A[Sequências FASTA] --> B[Alinhamento<br/>MAFFT/MUSCLE/pal2nal]
    B --> C[Filtragem<br/>trimAL]
    C --> D[Concatenação<br/>catfasta2phyml]
    D --> E[Inferência<br/>IQ-TREE 3 + ModelFinder]
    E --> F[Visualização<br/>iTOL]
    F --> G[Atribuição<br/>taxonômica]
```

Em projetos reais, essa mesma sequência se repete com variações de ferramenta, escala e organismo. **Os princípios são os mesmos** — o que muda são as flags e a interpretação contextual.

### Para se aprofundar

- **Métodos bayesianos (MrBayes, BEAST)** — incorporam incerteza nos parâmetros e permitem estimativas de tempo de divergência
- **Partições por códon** — modelar separadamente posições 1, 2 e 3 do códon pode melhorar a inferência
- **Filogenômica** — uso de centenas a milhares de genes para resolver relações profundas
- **Testes de topologia** — comparar hipóteses evolutivas alternativas (ex.: AU test, SH test)