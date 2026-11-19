# outputs/ — arquivos referenciados pelo tutorial de Nanopore

Esta pasta não é publicada como página, mas os arquivos dentro dela são linkados/embutidos
diretamente em `docs/cursos/nanopore/index.md` (relatórios de QC, imagens, tabelas de resultado).

Coloque os arquivos exatamente com os nomes abaixo (ou ajuste os links no `index.md` se preferir outros nomes).

## qc/
- `Scer1_NanoStats.txt`, `Scer2_NanoStats.txt`, `Scer1_NanoPlot-report.html`, `Scer2_NanoPlot-report.html` — NanoPlot dos dados brutos (Etapa 2)
- `exemplo_pico_DCS_lambda.png` — exemplo de pico de DCS/lambda (Etapa 2)
- `Scer1_q10_l500_NanoStats.txt`, `Scer2_q10_l500_NanoStats.txt`, `*_q10_l500_NanoPlot-report.html` — NanoPlot pós-filtragem (Etapa 5)

## kraken/
- `Scer1_kraken_report.txt`, `Scer2_kraken_report.txt` — relatório Kraken2 (Etapa 4, opcional)

## genomescope/
- `Scer1_Oxford_Nanopore_linear_plot.png`, `Scer2_Oxford_Nanopore_linear_plot.png` — GenomeScope2 (Etapa 5.2)

## bandage/
- `Scer1_Flye_bandage.png`, `Scer2_Flye_bandage.png` e os grafos `.gfa` — grafos de montagem (Etapa 6)

## quast/
- `report.html`, `icarus.html` — QUAST das montagens brutas (Etapa 9.1)

## merqury/
- `merqury_summary.tsv` — Merqury das montagens brutas (Etapa 9.2)
- `merqury_summary_polimento.tsv` — Merqury com k-mers Nanopore: bruta, purgada, Medaka e híbrida (Etapa 11.2)
- `merqury_illumina_summary.tsv` — Merqury com k-mers Illumina, só Scer1 (Etapa 11.3)

## busco/
- `short_summary.specific.saccharomycetaceae_odb12.<montagem>.txt` — BUSCO das montagens brutas (Etapa 9.3)
- `busco_figure.png` — gráfico comparativo (Etapa 9.3)

## purge_dups/
- `<montagem>.cobertura.png` — histogramas de cobertura do purge_dups, Scer1 (Etapa 9.5)

## quast-polimento/
- `report.html`, `icarus.html` — QUAST das 16 versões (bruta, purgada, Medaka, híbrida) (Etapa 11.1)

## dotplot/
- `dotplot_Scer1_Flye_{contigs,ragtag,longstitch}.png` — dotplots contra a S288C (Etapa 12.5)

## blobplot/
- `<montagem>.blob.png` — blob plots do BlobToolKit (Etapa 13.3)

## telomeros/
- `tidk_Scer1_hifiasm.svg` — `tidk plot` da Scer1 Hifiasm final (Etapa 14.3)

## snailplot/
- `<montagem>.snail.png` — snail plots do BlobToolKit (Etapa 14.4)

---

Os FASTAs de montagem (grandes) **não** vão aqui — eles ficam em
`data/exemplos/nanopore/assemblies/` na raiz do repositório (ver README lá).
