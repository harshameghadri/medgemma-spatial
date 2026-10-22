# Knowledge snapshots

`omnipath_lr.tsv.gz` (ligand-receptor pairs with signs) and `omnipath_roles.tsv.gz`
(ligand / receptor role per gene) are derived from OmniPath by
`scripts/fetch_omnipath.py`; `omnipath_lr.json` records the exact queries and the
fetch date (2026-10-22).

- Source: OmniPath, https://omnipathdb.org (`interactions` with datasets
  `omnipath,ligrecextra`; `intercell` categories `ligand,receptor`).
- Only resources that OmniPath marks as allowing commercial use were requested
  (`license=commercial`). The licences of the underlying resources still apply to
  the interactions; see https://omnipathdb.org/info.
- Cite: Türei D, Korcsmáros T, Saez-Rodriguez J. OmniPath: guidelines and gateway
  for literature-curated signaling pathway resources. *Nat Methods* 13, 966–967 (2016).
  Türei D et al. Integrated intra- and intercellular signaling knowledge for
  multicellular omics analysis. *Mol Syst Biol* 17, e9923 (2021).
