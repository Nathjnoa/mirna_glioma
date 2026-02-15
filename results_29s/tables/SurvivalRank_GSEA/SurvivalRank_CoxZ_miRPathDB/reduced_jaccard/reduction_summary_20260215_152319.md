# Pathway Redundancy Reduction Summary

**Method**: Jaccard similarity + hierarchical clustering (Enrichment Map)
**Reference**: Merico D et al. (2010) PLoS ONE 5(11):e13984

- **Run tag**: SurvivalRank_CoxZ_miRPathDB
- **Linkage method**: average (UPGMA)
- **Date**: 2026-02-15 15:40

## Database-specific parameters:

**KEGG pathways**:
  - Jaccard threshold: 0.25 (permissive, due to sparse overlap)
  - Min cluster size: 1 (retain singletons)

**Reactome pathways**:
  - Jaccard threshold: 0.5 (standard)
  - Min cluster size: 2

**Rationale**: KEGG pathways exhibit lower inherent gene overlap and smaller
pathway count compared to Reactome, warranting a more permissive threshold to
avoid over-reduction while maintaining biological interpretability.

## Overall reduction:
- **Before reduction**: 2719
- **After reduction**: 2132
- **Removed**: 587 redundant pathways (21.6%)

## Reduction per database/direction

| Database | Direction | Method | Input | Output | Reduction % | Clusters | Median Jaccard | Threshold | Min Size |
|----------|-----------|--------|-------|--------|-------------|----------|----------------|-----------|----------|
| KEGG | enriched | Jaccard_hclust | 24 | 14 | 41.7% | 14 | 0.042 | 0.25 | 1 |
| KEGG | depleted | Jaccard_hclust | 105 | 36 | 65.7% | 36 | 0.132 | 0.25 | 1 |
| Reactome | enriched | Jaccard_hclust | 129 | 19 | 85.3% | 19 | 0.040 | 0.5 | 2 |
| Reactome | depleted | Jaccard_hclust | 518 | 93 | 82.0% | 93 | 0.048 | 0.5 | 2 |
| GO_BP | enriched | rrvgo_semantic | 183 | 17 | 90.7% | - | - | 0.9 | - |
| GO_BP | depleted | rrvgo_semantic | 1276 | 31 | 97.6% | - | - | 0.9 | - |

## Pathways per category (after reduction)

| Database | Enriched | Depleted |
|----|----:|----:|
| GO_BP | 219 | 1505 |
| KEGG | 14 | 35 |
| Reactome | 69 | 290 |
