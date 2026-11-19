# igv/ — arquivos para visualizar as junções do chrV no IGV (Etapa 15)

Recortes da montagem Scer1 Hifiasm final e das reads Nanopore alinhadas a ela (BAM da Etapa 13.2):

- `hifiasm_chrV.fasta` (+ `.fai`) — contigs `ptg000007l_1` e `ptg000023l_1`
- `hifiasm_juncoes.bam` (+ `.bai`) — reads de `ptg000007l_1:190000-215000` (controle), `ptg000007l_1:430000-470000` (junção chrV|chrXVI) e `ptg000023l_1:110000-145000` (junção chrV|chrXIV)

Para abrir em https://igv.org/app: *Genome → Local File* (fasta + fai) e *Tracks → Local File* (bam + bai).
