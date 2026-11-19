# assemblies/ — FASTAs do tutorial de Nanopore

FASTAs referenciados em `docs/cursos/nanopore/index.md` (links diretos para o GitHub, fora do build do
mkdocs por serem arquivos grandes). Nomes esperados:

**Montagens brutas**
- `Scer1_Flye_assembly.fasta`, `Scer2_Flye_assembly.fasta`
- `Scer1_Hifiasm_assembly.fa`, `Scer2_Hifiasm_assembly.fa`
- `Scer1_NextDenovo_assembly.fasta`

**Polidas com Medaka** (só a parte nuclear, sem a mitocôndria; na Scer1, também depois do purge_dups)
- `Scer1_Flye_medaka_consensus.fasta`, `Scer2_Flye_medaka_consensus.fasta`
- `Scer1_Hifiasm_medaka_consensus.fasta`, `Scer2_Hifiasm_medaka_consensus.fasta`
- `Scer1_NextDenovo_medaka_consensus.fasta`

**Scaffolded** (a partir das versões de `9.Montagens-polidas/`: Scer1 híbrida, Scer2 Medaka)
- `Scer1_Flye_ragtag_scaffold.fasta`, `Scer2_Flye_ragtag_scaffold.fasta`, `Scer1_Hifiasm_ragtag_scaffold.fasta`
- `Scer1_Flye_longstitch_scaffold.fasta`, `Scer2_Flye_longstitch_scaffold.fasta`

Se preferir outros nomes, ajuste os links correspondentes em `docs/cursos/nanopore/index.md`.
