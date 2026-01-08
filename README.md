# mirna_glioma

Analisis de expresion de miRNA en glioma con controles de calidad, expresion diferencial,
supervivencia y enriquecimiento funcional.

## Estructura del repo
- `config/`: parametros de analisis (por ejemplo, `config/de_specs.csv`).
- `docs/`: documentacion metodologica y notas de ejecucion.
- `env/`: definicion del entorno reproducible.
- `scripts/`: pipeline en R (QC, DE, heatmaps, supervivencia, GSEA).

Carpetas generadas localmente (ignoradas por Git):
- `data/`: entradas locales y datos intermedios/processed.
- `results/`: tablas y figuras de salida.
- `logs/`: logs con timestamp por ejecucion.

## Inicio rapido
1) Crear el entorno:
```bash
conda env create -f env/environment.yml
```
2) Revisar y ajustar `config/de_specs.csv` segun las comparaciones deseadas.
3) Ejecutar scripts principales en orden (ver `scripts/`).

## Documentacion
- Metodos y descripcion del pipeline: `docs/METHODS_scripts.md`.

## Reproducibilidad
- Cada script genera logs con timestamp en `logs/`.
- Las salidas se escriben en `results/` y `data/processed/`.

## Licencia
- {PENDIENTE}
