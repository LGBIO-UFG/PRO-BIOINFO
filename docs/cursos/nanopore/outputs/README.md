# outputs/ — arquivos referenciados pelo tutorial de Nanopore

Esta pasta não é publicada como página, mas os arquivos dentro dela são linkados/embutidos
diretamente em `docs/cursos/nanopore/index.md` (relatórios de QC, imagens, tabelas de resultado).

Coloque os arquivos exatamente com os nomes abaixo (ou ajuste os links no `index.md` se preferir outros nomes).

## qc/
- `Scer1_NanoStats.txt`, `Scer2_NanoStats.txt` — estatísticas NanoPlot dos dados brutos (Etapa 2)
- `Scer1_NanoPlot-report.html`, `Scer2_NanoPlot-report.html` — relatório interativo NanoPlot dos dados brutos
- `Scer1_q10_l500_NanoStats.txt`, `Scer2_q10_l500_NanoStats.txt` — estatísticas NanoPlot pós-filtragem (Etapa 5)
- `Scer1_q10_l500_NanoPlot-report.html`, `Scer2_q10_l500_NanoPlot-report.html` — relatório NanoPlot pós-filtragem

## kraken/
- `Scer1_kraken_report.txt`, `Scer2_kraken_report.txt` — relatório taxonômico Kraken2 (Etapa 4, opcional)

## genomescope/
- `Scer1_Oxford_Nanopore_linear_plot.png`, `Scer2_Oxford_Nanopore_linear_plot.png` — gráfico linear GenomeScope2 (Etapa 5.2, embutido na página)

## quast/
- `report.html`, `icarus.html` — QUAST das montagens brutas (Etapa 9.1)

## quast-medaka/
- `report.html`, `icarus.html` — QUAST pós-polimento com Medaka (Etapa 11.1)

## merqury/
- `merqury_summary.tsv` — resumo Merqury das montagens brutas (Etapa 9.2)
- `merqury_summary_medaka.tsv` — resumo Merqury pós-polimento (Etapa 11.2)

## busco/
- `short_summary.specific.saccharomycetaceae_odb12.Scer1_q10_l500_Flye.txt`
- `short_summary.specific.saccharomycetaceae_odb12.Scer1_q10_l500_hifiasm.txt`
- `short_summary.specific.saccharomycetaceae_odb12.Scer1_q10_l500_nextdenovo.txt`
- `short_summary.specific.saccharomycetaceae_odb12.Scer2_q10_l500_Flye.txt`
- `short_summary.specific.saccharomycetaceae_odb12.Scer2_q10_l500_hifiasm.txt`
- `busco_figure.png` — gráfico comparativo gerado pelo `generate_plot.py` (Etapa 9.3, embutido na página)

## blobtools/
- `Scer1_q10_l500_Flye.blobplot.genus.png`, `Scer2_q10_l500_Flye.blobplot.genus.png` — blob plots (Etapa 14, embutidos na página)
- `Scer1_q10_l500_Flye.blobDB.table.txt`, `Scer2_q10_l500_Flye.blobDB.table.txt` — tabela BlobTools por contig

---

Os FASTAs de montagem (grandes) **não** vão aqui — eles ficam em
`data/exemplos/nanopore/assemblies/` na raiz do repositório (ver README lá).
