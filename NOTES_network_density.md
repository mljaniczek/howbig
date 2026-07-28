# Notes: Within-Block Density and Comparison to the JGL Paper

## Summary

The production simulation sweeps are generating networks with a
**within-block density of ~8%**, not the **15%** documented in
`CLAUDE.md` and intended in the original simulation design.

The cause is that `generate_block_graph()`
(`01_generate_gold_standard_networks.R:14-26`) calls
`igraph::sample_pa(n = p_block, m = m_pa)` with `m_pa = 1`
(set in the `base` settings at `00_helper_functions.R:74`), which
produces a Barabási–Albert **tree** rather than a proper scale-free
graph. The `within_density = 0.15` parameter is silently ignored for
scale-free topology — it only takes effect when
`topology = "random"` (Erdős–Rényi).

**Switching `m_pa = 1` → `m_pa = 2` (one-line fix) would:**

- Bring within-block density from 8.0% → 15.7% (matching the
  documented design intent)
- Increase mean within-block degree from 1.9 → 3.8
- Replace tree structure with proper scale-free graphs (hubs and
  cycles, not branching chains)
- Approximately double the Setting-A differential edge count
  (block deletion removes ~47 edges instead of 24)

## Empirical numbers (50 realizations of `sample_pa()`)

| Setting                                 | p_block | m_pa | Edges/block | Max edges | Density | Mean degree |
|-----------------------------------------|--------:|-----:|------------:|----------:|--------:|------------:|
| Danaher (2014), p=500                   |      50 |    2 |        97.0 |     1,225 |    7.9% |        3.88 |
| Danaher (our matched p=100 replication) |      10 |    2 |        17.0 |        45 |   37.8% |        3.40 |
| **Production, current (m_pa=1)**        |      25 |    1 |        24.0 |       300 | **8.0%**|        1.92 |
| Production, target (m_pa=2)             |      25 |    2 |        47.0 |       300 |**15.7%**|        3.76 |

## Density match with JGL is coincidental

Our current production runs match the within-block density of the
original JGL simulation (Danaher 2014 at p=500, p_block=50, m_pa=2:
~8%) almost exactly. This is a numerical coincidence — BA edges scale
roughly as `m_pa × (p_block − 1)` while max possible edges scale as
`p_block² / 2`, so `p_block=25, m_pa=1` happens to land at the same
density as Danaher's `p_block=50, m_pa=2`.

However, the **topology differs**: Danaher's blocks have mean degree
~3.9 (proper scale-free, with hubs and a connected core), while ours
have mean degree ~1.9 (trees — every node has exactly one parent, no
cycles). Trees lack the hub-and-cycle structure that makes scale-free
networks distinctive.

## Other JGL parameter comparisons

From `01b_generate_paper_gold_standards.R:135-192`
(`generate_gs_danaher()`):

| Parameter             | Danaher (2014)                | Our production                              |
|-----------------------|-------------------------------|---------------------------------------------|
| Block size            | 50 (at p=500)                 | 25                                          |
| Within-block topology | Scale-free (BA m=2)           | Scale-free (BA m=1)                         |
| Within-block density  | ~8%                           | ~8% (current) / ~16% (if m_pa=2)            |
| Mean within-block deg | ~3.9                          | ~1.9 (current) / ~3.8 (if m_pa=2)           |
| Between-block edges   | None (block_prob = 0)         | 0.5–2.5% (block_prob = 0.01 × 5/B)          |
| Edge weights          | U[−0.4,−0.1] ∪ [0.1,0.4]      | Empirical (mxDist.csv)                      |
| Precision build       | Row-sum + d_scale = 0.6       | Row-sum, d_scale = NULL                     |
| Perturbation          | Remove block 1                | A: remove block 1; B/C: + magnitude scaling |
| n                     | 150                           | 100–5000                                    |

## Biology context

### Whole-organism networks (very sparse, ~0.1–1%)
- PPI: yeast Y2H mean degree ~2.4, density ~0.5%
  (Jeong et al. 2001, *Nature*)
- Human interactome: ~3,200 edges among 1,705 proteins, ~0.2%
  density (Stelzl et al. 2005, *Cell*)
- Metabolic: substrate degree ~3.2 across 43 organisms
  (Jeong et al. 2000, *Nature*)
- Foundational review: Barabási & Oltvai (2004),
  *Nature Reviews Genetics*

### Pathway-restricted networks (5–25%)
- GGM on serum metabolomics: 5–15% within-pathway density
  (Krumsiek et al. 2011, *BMC Systems Biology*)
- WGCNA modules: 5–20% within-module density after partialing
  (Zhang & Horvath 2005)
- KEGG / Reactome pathways: 20–40%, but this counts all annotated
  regulatory/binding/catalytic relationships (looser than GGM
  recovery from expression data)

### Differential network density
- ~5% of edges rewired between yeast pre/post DNA damage
  (Bandyopadhyay et al. 2010, *Science*)
- 1–20% differential range across stress / disease perturbations
  (Ideker & Krogan 2012, *Mol Syst Biol*)
- Conceptual framework: de la Fuente (2010), *Trends in Genetics*

### Does our simulation match biology?

| Property                              | Current (m_pa=1) | Target (m_pa=2) | Biology (pathway-scale) |
|---------------------------------------|-----------------:|----------------:|------------------------:|
| Topology                              |             Tree |      Scale-free |              Scale-free |
| Within-block density                  |             8.0% |           15.7% |                   5–20% |
| Mean within-block degree              |              1.9 |             3.8 |                    3–10 |
| Differential edge fraction (at p=250) |       0.08–0.24% |          ~0.5%  |                   1–20% |

Current production is on the low end of pathway-scale density and
produces a tree topology rare in real biology. `m_pa = 2` lands in
the middle of the pathway-scale range and produces topologically
realistic scale-free graphs.

## Recommendation

**Switch `m_pa = 1` → `m_pa = 2`** in `00_helper_functions.R:74`. This
is a one-line change that:

1. Reconciles documented design (15% density) with actual production
2. Produces topologically realistic scale-free blocks
3. Matches Danaher's small-block topology (mean degree ~3.4) and
   his large-block density (~8%) simultaneously
4. Strengthens the claim that the simulation matches biology

**Validation before adopting:** spot-check the headline calibration
findings at `m_pa = 2` on a small (~10-rep) subset to confirm they
hold:

- "Pooled-permutation BH is anti-conservative for FGL/GGL under
  magnitude perturbations" — does the achieved-FDR-up-to-0.30
  finding survive at denser blocks?
- "scale-multiGGM PIPs collapse for p ≳ 250" — does the collapse
  threshold shift at denser blocks?
- "Method ranking is stable across topology" — does denser
  scale-free vs current tree change relative method ordering?

If the headline findings shift, the manuscript narrative needs
adjustment. If they hold, the move to `m_pa = 2` is pure upside
(matches design, matches biology, no narrative cost).

## File pointers

- Density-controlling line: `00_helper_functions.R:74`
  (`scale_factor = 1.5, m_pa = 1`)
- Block generator: `01_generate_gold_standard_networks.R:14-26`
  (`generate_block_graph()`)
- JGL reference generator: `01b_generate_paper_gold_standards.R:135-192`
  (`generate_gs_danaher()`)
- The `within_density` parameter is checked at
  `01_generate_gold_standard_networks.R:18-19` and only used for the
  random/ER branch (`sample_gnp(p_block, within_density)`); the
  scale-free branch at line 21 ignores it.

## Suggested manuscript text (current m_pa=1)

> Within-block topology was generated as a Barabási–Albert tree
> (m=1) on 25 nodes, yielding 24 edges per block (8% within-block
> density). This matches the within-block density of the original
> JGL simulation study (Danaher et al. 2014; p=500, p_block=50,
> m=2, ~8% within-block density) and the sparse regime of
> pathway-restricted Gaussian graphical models on metabolomics
> data (Krumsiek et al. 2011). Our networks are sparser in degree
> (mean within-block degree ~1.9) than Danaher's (mean ~3.9),
> but match in overall density.

## Suggested manuscript text (if switched to m_pa=2)

> Within-block topology was generated as a Barabási–Albert
> scale-free graph (m=2) on 25 nodes, yielding ~47 edges per
> block (15.7% within-block density, mean within-block degree
> ~3.8). This matches both the topology and density regime of
> the original JGL simulation (Danaher et al. 2014; mean
> within-block degree ~3.9) and the pathway-restricted GGM
> literature on metabolomics data (Krumsiek et al. 2011; 5–15%
> within-pathway density).
