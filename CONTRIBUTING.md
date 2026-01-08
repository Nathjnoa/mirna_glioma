# Contribucion

Gracias por querer contribuir a `mirna_glioma`. Este repo es un pipeline de bioinformatica
en R para analisis de expresion de miRNA en glioma.

## Reglas basicas
- No subir datos sensibles ni muestras identificables.
- No commitear `data/`, `results/` ni `logs/` (estan ignorados por defecto).
- Mantener cambios pequenos y con mensajes claros.

## Entorno y ejecucion
1) Crear el entorno:
```bash
conda env create -f env/environment.yml
```
2) Revisar y ajustar `config/de_specs.csv`.
3) Ejecutar scripts desde la raiz del repo, por ejemplo:
```bash
Rscript scripts/01_inspeccion_counts.R
```

## Estilo y buenas practicas
- Evitar cambiar la logica cientifica sin documentar el motivo.
- Si agregas un parametro nuevo, documentalo en `README.md` y/o `docs/METHODS_scripts.md`.
- Mantener salidas reproducibles (logs con timestamp ya estan en los scripts).

## Flujo de Git recomendado
```bash
git status -sb
git add <archivo>
git commit -m "tipo: descripcion corta"
git push
```

Ejemplos de mensajes:
- `feat: add new QC plot`
- `fix: handle missing metadata`
- `docs: update METHODS`

## Issues y cambios grandes
- Abre un issue para cambios mayores (nuevos modulos, nuevos datasets).
- Explica el objetivo, los pasos y el impacto esperado.
