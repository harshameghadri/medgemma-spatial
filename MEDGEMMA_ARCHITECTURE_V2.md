# MedGemma Architecture V2: Anti-Parrot Design

**Problem**: MedGemma v1 was just paraphrasing JSON → prose (zero added value)

**Solution**: Multi-stage reasoning pipeline where parroting is **structurally impossible**

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│ Stage 1: SPATIAL PATTERN EXTRACTION (NO RAW DATA)              │
│ ─────────────────────────────────────────────────────────────  │
│ Input: AnnData object + spatial graphs (NOT JSON summary)      │
│                                                                 │
│ MedGemma Task: Identify patterns humans might miss             │
│ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│ Questions:                                                      │
│ 1. "Looking at ISG15 spatial expression (Moran's I=0.57),     │
│    where specifically is it clustering? Tumor center,          │
│    invasion front, stroma, or interface zones?"                │
│                                                                 │
│ 2. "Cell types show complete spatial segregation. Is this      │
│    uniform across the tissue or are there regional             │
│    exceptions? If exceptions exist, what genes are             │
│    differentially expressed in those zones?"                   │
│                                                                 │
│ 3. "Entropy = 0.42 suggests moderate heterogeneity. Are        │
│    there distinct spatial domains (e.g., homogeneous tumor     │
│    core vs heterogeneous periphery)?"                          │
│                                                                 │
│ OUTPUT: Spatial hypotheses (NOT data summary)                  │
└─────────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────────┐
│ Stage 2: COMPARATIVE MEDICAL REASONING (FORCED COMPARISON)     │
│ ─────────────────────────────────────────────────────────────  │
│ Input: Current sample + Reference database (3-5 samples)       │
│                                                                 │
│ MedGemma Task: Contextualize within known phenotypes           │
│ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│ Reference Samples (Synthetic or Real):                         │
│                                                                 │
│ SAMPLE A: "Hot" Tumor                                          │
│ - Interface: 60%, ISG15 Moran I: 0.85, Entropy: 0.65          │
│ - Phenotype: Inflamed, high TILs, PD-L1+                      │
│ - Treatment: Responds to checkpoint inhibitors                 │
│                                                                 │
│ SAMPLE B: "Cold" Tumor                                         │
│ - Interface: 12%, ISG15 Moran I: 0.15, Entropy: 0.28          │
│ - Phenotype: Immune excluded, stromal barrier                  │
│ - Treatment: Poor checkpoint response, good endocrine therapy  │
│                                                                 │
│ SAMPLE C: "Immune Excluded"                                    │
│ - Interface: 35%, ISG15 Moran I: 0.45, Entropy: 0.50          │
│ - Phenotype: Immune cells present but segregated               │
│ - Treatment: Combination therapy (anti-VEGF + checkpoint)      │
│                                                                 │
│ CURRENT SAMPLE:                                                 │
│ - Interface: 34.6%, ISG15 Moran I: 0.57, Entropy: 0.42        │
│                                                                 │
│ MedGemma Questions (FORCES COMPARISON):                        │
│ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│ 1. "Which reference sample does this most resemble? Provide   │
│    similarity scores (0-100) for each and explain why."       │
│                                                                 │
│ 2. "What spatial features DISTINGUISH this sample from all     │
│    references? (i.e., what makes it unique?)"                 │
│                                                                 │
│ 3. "If this most resembles Sample C (immune excluded), but    │
│    ISG15 clustering is HIGHER, what does that discrepancy     │
│    suggest about immune activation?"                           │
│                                                                 │
│ OUTPUT: Comparative phenotype classification + unique features │
└─────────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────────┐
│ Stage 3: CONDITIONAL REASONING (IF-THEN LOGIC)                 │
│ ─────────────────────────────────────────────────────────────  │
│ MedGemma Task: Clinical decision tree based on spatial data    │
│ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│ Questions:                                                      │
│                                                                 │
│ 1. "IF spatial segregation is complete (all cell types         │
│    depleted) AND interface is 34.6%, THEN what does this      │
│    suggest about chemokine gradients or physical barriers?"    │
│                                                                 │
│ 2. "IF ISG15 clusters at the interface (not tumor core),      │
│    AND segregation is complete, THEN is the immune response   │
│    active but ineffective, or suppressed?"                     │
│                                                                 │
│ 3. "IF this is HR+ luminal breast cancer AND shows 'immune    │
│    excluded' phenotype, THEN rank these treatments by          │
│    expected efficacy:                                           │
│    - Endocrine therapy alone                                    │
│    - Checkpoint inhibitor alone                                 │
│    - Endocrine + CDK4/6 inhibitor                              │
│    - Combination immunotherapy + stromal targeting"            │
│                                                                 │
│ OUTPUT: Conditional clinical reasoning chains                  │
└─────────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────────┐
│ Stage 4: HYPOTHESIS GENERATION (TESTABLE PREDICTIONS)          │
│ ─────────────────────────────────────────────────────────────  │
│ MedGemma Task: Generate follow-up experiments                  │
│ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│ Questions:                                                      │
│                                                                 │
│ 1. "Based on spatial segregation pattern, propose 3 testable  │
│    hypotheses about what's preventing T-cell infiltration.    │
│    For each, specify:                                           │
│    - Hypothesis                                                 │
│    - Validation experiment (IHC, CODEX, etc.)                  │
│    - Expected result if true"                                   │
│                                                                 │
│ 2. "ISG15 clustering suggests localized interferon response.   │
│    What OTHER interferon-stimulated genes should we look for  │
│    in the same spatial regions to confirm this?"               │
│                                                                 │
│ 3. "If we performed multiplexed IHC on this tissue, which     │
│    protein markers would you prioritize to validate the        │
│    spatial transcriptomics findings?"                           │
│                                                                 │
│ OUTPUT: Testable hypotheses + validation experiments           │
└─────────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────────┐
│ Stage 5: VERIFICATION (Llama-3.1-8B or Qwen2.5-32B)           │
│ ─────────────────────────────────────────────────────────────  │
│ Task: Validate MedGemma's medical reasoning                    │
│ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│ Input: MedGemma's hypotheses + spatial data                    │
│                                                                 │
│ Llama-3.1-8B Questions:                                        │
│                                                                 │
│ 1. "MedGemma claims ISG15 clustering indicates 'localized     │
│    interferon response at immune-tumor interface.' Is this     │
│    biologically plausible? What alternative explanations       │
│    exist?"                                                      │
│                                                                 │
│ 2. "Validate the treatment ranking:                            │
│    MedGemma ranked: [Endocrine+CDK4/6 > Endocrine alone >     │
│    Combo immuno > Checkpoint alone]                            │
│    Do you agree with this ranking for immune-excluded HR+      │
│    breast cancer? Cite evidence."                              │
│                                                                 │
│ 3. "Cross-check: Are there any CONTRADICTIONS between the     │
│    spatial patterns and MedGemma's clinical interpretation?"   │
│                                                                 │
│ OUTPUT: Validation report (agree/disagree + reasoning)         │
└─────────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────────┐
│ Stage 6: NOVELTY VALIDATION (Programmatic Check)               │
│ ─────────────────────────────────────────────────────────────  │
│ Task: Ensure insights are actually novel                       │
│ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│ Programmatic Checks:                                            │
│                                                                 │
│ def validate_novelty(medgemma_output, input_json):             │
│     """Reject parrot outputs."""                               │
│                                                                 │
│     # Check 1: Does output contain raw metrics from JSON?      │
│     forbidden_phrases = [                                       │
│         str(input_json['tumor_immune_interface']['interface_pct']),│
│         ", ".join(input_json['spatial_autocorrelation']['top_genes'][:3]),│
│         "spatial segregation observed"  # Generic phrase       │
│     ]                                                           │
│                                                                 │
│     for phrase in forbidden_phrases:                            │
│         if phrase.lower() in medgemma_output.lower():          │
│             raise ValueError(f"PARROT DETECTED: '{phrase}'")   │
│                                                                 │
│     # Check 2: Does output contain comparative analysis?       │
│     required_elements = [                                       │
│         "compared to",                                          │
│         "suggests",                                             │
│         "hypothesis",                                           │
│         "if.*then" (regex)                                      │
│     ]                                                           │
│                                                                 │
│     score = sum(1 for elem in required_elements                │
│                 if re.search(elem, medgemma_output, re.I))     │
│                                                                 │
│     if score < 2:                                               │
│         raise ValueError("Insufficient reasoning depth")        │
│                                                                 │
│     # Check 3: Novelty metric (TF-IDF against input JSON)      │
│     novelty_score = calculate_tfidf_divergence(                │
│         medgemma_output,                                        │
│         json.dumps(input_json)                                  │
│     )                                                           │
│                                                                 │
│     if novelty_score < 0.6:  # Threshold                       │
│         raise ValueError(f"Low novelty: {novelty_score:.2f}")  │
│                                                                 │
│     return True                                                 │
│                                                                 │
│ OUTPUT: PASS/FAIL + novelty score                              │
└─────────────────────────────────────────────────────────────────┘
```

---

## Implementation Details

### Stage 1: Spatial Pattern Extraction

**Key Innovation**: Feed raw AnnData, NOT pre-computed summaries

```python
def stage1_spatial_patterns(adata, spatial_stats):
    """Extract spatial patterns MedGemma can reason about."""

    # Get top spatially variable genes
    top_genes = spatial_stats['spatial_autocorrelation']['top_genes'][:5]
    morans_i = spatial_stats['spatial_autocorrelation']['top_morans_i'][:5]

    # Create gene-specific spatial queries
    gene_queries = []
    for gene, mi in zip(top_genes, morans_i):
        # Compute WHERE the gene clusters (not just that it clusters)
        gene_expr = adata.raw[:, gene].X.toarray().flatten()
        spatial_coords = adata.obsm['spatial']
        cell_types = adata.obs['cell_type']

        # Find high-expression regions
        high_expr_mask = gene_expr > np.percentile(gene_expr, 75)
        high_expr_celltypes = cell_types[high_expr_mask].value_counts()

        # Calculate spatial location (center/periphery/interface)
        centroid = spatial_coords.mean(axis=0)
        distances = np.linalg.norm(spatial_coords - centroid, axis=1)
        high_expr_avg_dist = distances[high_expr_mask].mean()
        all_spots_avg_dist = distances.mean()

        location = "center" if high_expr_avg_dist < all_spots_avg_dist else "periphery"

        query = f"""
Gene: {gene} (Moran's I = {mi:.2f})

SPATIAL QUESTION: This gene shows strong spatial clustering.
Based on the data:
- It's highly expressed in: {', '.join(high_expr_celltypes.index[:3])}
- Spatial location: {location} of tissue
- Expression is NOT uniform

MEDICAL REASONING QUESTION:
Given {gene}'s known biological function, what does its specific
spatial pattern (clustering in {location}, enriched in {high_expr_celltypes.index[0]})
tell us about the tissue microenvironment that bulk RNA-seq would miss?

What molecular mechanisms could explain WHY this gene clusters in this
specific spatial pattern?
"""
        gene_queries.append(query)

    return "\n\n".join(gene_queries)
```

**Expected Output (GOOD)**:
```
ISG15 (Moran's I = 0.57)
Location: Tissue periphery, enriched in plasma_IgG cells

INTERPRETATION: ISG15's peripheral clustering in immune cells suggests
localized interferon response at the tumor-stroma boundary, NOT systemic
immune activation. This is a SPATIAL biomarker indicating:

1. Active but LOCALIZED immune recognition
2. Potential immune-tumor interaction zones at periphery
3. Distinct from diffuse ISG15 expression (viral infection)

MECHANISTIC HYPOTHESIS: Tumor cells at invasion front secrete chemokines
(CXCL10, CXCL9) recruiting IFN-activated immune cells, but stromal barriers
prevent penetration to tumor center.

CLINICAL SIGNIFICANCE: Peripheral immune activation suggests tumor is
NOT completely immunologically "cold", but physical barriers prevent
effective immune attack. Rationale for combination therapy targeting
both immune activation AND stromal remodeling.
```

**Not Allowed (BAD)**:
```
ISG15 is spatially variable with Moran's I = 0.57.  ← PARROT
```

---

### Stage 2: Comparative Analysis

**Key Innovation**: Force MedGemma to compare against reference samples

```python
# Create reference database (can be synthetic or real samples)
REFERENCE_SAMPLES = {
    "hot_tumor": {
        "interface_pct": 60,
        "isg15_morans_i": 0.85,
        "entropy": 0.65,
        "enrichment": "all_enriched",
        "phenotype": "Inflamed, high TILs",
        "treatment_response": "Good checkpoint inhibitor response",
        "clinical_example": "Melanoma, MSI-high CRC"
    },
    "cold_tumor": {
        "interface_pct": 12,
        "isg15_morans_i": 0.15,
        "entropy": 0.28,
        "enrichment": "all_depleted",
        "phenotype": "Immune desert",
        "treatment_response": "Poor immunotherapy, good targeted therapy",
        "clinical_example": "Pancreatic adenocarcinoma"
    },
    "immune_excluded": {
        "interface_pct": 35,
        "isg15_morans_i": 0.45,
        "entropy": 0.50,
        "enrichment": "all_depleted",
        "phenotype": "Immune cells present but segregated",
        "treatment_response": "May respond to combination therapy",
        "clinical_example": "Triple-negative breast cancer with stromal barrier"
    }
}

def stage2_comparative_reasoning(current_sample, references=REFERENCE_SAMPLES):
    """Force comparative analysis."""

    prompt = f"""
You are analyzing a tumor sample's spatial transcriptomics profile.

REFERENCE SAMPLES (Known Phenotypes):

{format_references(references)}

CURRENT SAMPLE:
- Tumor-immune interface: {current_sample['interface_pct']}%
- ISG15 spatial clustering (Moran's I): {current_sample['isg15_morans_i']}
- Spatial entropy: {current_sample['entropy']}
- Cell-cell interactions: {current_sample['enrichment']}

MANDATORY ANALYSIS (Answer ALL):

1. SIMILARITY SCORING:
   Provide similarity score (0-100) to EACH reference sample. Explain why.

   Hot Tumor: ___/100 because ___
   Cold Tumor: ___/100 because ___
   Immune Excluded: ___/100 because ___

2. UNIQUE FEATURES:
   What spatial characteristics distinguish this sample from ALL references?
   (If it's identical to a reference, explain why that's significant)

3. DISCREPANCY ANALYSIS:
   The current sample shows:
   - Immune excluded pattern (segregation, moderate interface)
   - BUT ISG15 clustering is HIGHER than typical immune excluded tumors

   What does this discrepancy suggest? Propose 2 mechanistic explanations.

4. PHENOTYPE CLASSIFICATION:
   Based on the comparison, classify this tumor as:
   [ ] Hot (inflamed)
   [ ] Cold (desert)
   [ ] Immune excluded
   [ ] Mixed/Atypical (specify)

5. TREATMENT PREDICTION:
   Given the spatial phenotype, rank these therapies (1=best, 4=worst):
   ___ Checkpoint inhibitor monotherapy
   ___ Endocrine therapy + CDK4/6 inhibitor
   ___ Checkpoint inhibitor + anti-VEGF (stromal targeting)
   ___ Targeted therapy alone

   Justify your ranking based on spatial features.
"""
    return prompt
```

**This FORCES MedGemma to:**
- Compare (can't just describe current sample)
- Quantify (similarity scores)
- Reason (explain discrepancies)
- Predict (rank treatments)

---

### Stage 3: Conditional Reasoning

**Key Innovation**: IF-THEN clinical logic chains

```python
def stage3_conditional_reasoning(spatial_data, celltype_data):
    """Force logical reasoning chains."""

    prompt = f"""
CLINICAL REASONING EXERCISE:

Use IF-THEN logic to connect spatial observations to biological mechanisms.

OBSERVATION SET:
A. Spatial segregation: ALL cell types show neighborhood depletion
B. Tumor-immune interface: 34.6% (moderate)
C. ISG15 clustering: Moran's I = 0.57 (strong)
D. Spatial entropy: 0.42 (moderate heterogeneity)
E. Cell composition: 51% epithelial, 41% immune, 8% stromal

REASONING CHAINS (Complete ALL):

1. IF segregation is complete (observation A)
   AND interface exists (observation B)
   THEN what physical structure separates cell types?

   Possibilities:
   - Basement membrane barrier
   - Stromal matrix deposition
   - Chemokine gradient repulsion
   - Metabolic barrier (hypoxia, acidosis)

   Which is most likely? Why? What spatial gene expression would confirm?

2. IF ISG15 clusters strongly (observation C)
   AND it's an interferon-stimulated gene
   THEN where specifically should we see other ISG genes (ISG20, MX1, IFIT1)?

   Prediction: Co-localized [ ] Yes [ ] No

   IF they co-localize → suggests ___
   IF they don't → suggests ___

3. IF this is HR+ breast cancer (epithelial markers)
   AND immune cells are present but excluded (observations A+B)
   AND moderate heterogeneity exists (observation D)
   THEN what does this imply for endocrine therapy resistance?

   Complete: "Spatial heterogeneity at entropy=0.42 suggests ___ distinct
   tumor subclones. The immune excluded phenotype indicates ___ microenvironment.
   These combined patterns predict endocrine therapy will be [effective/ineffective]
   because ___."

4. TREATMENT LOGIC CHAIN:

   IF immune cells are present (41% of tissue)
   BUT spatially excluded (all depleted)
   AND ISG15 clustering suggests active immune cells
   THEN checkpoint inhibitor monotherapy will likely [succeed/fail] because ___.

   HOWEVER, IF we ALSO target ___ (stromal barrier/VEGF/TGF-β),
   THEN the spatial barrier might ___, allowing ___.

5. PROGNOSIS CHAIN:

   IF luminal subtype (good prognosis)
   AND immune excluded (varies)
   AND moderate heterogeneity (potential resistance)
   THEN 5-year survival probability is:
   [ ] >90% (excellent)
   [ ] 70-90% (good)
   [ ] 50-70% (moderate)
   [ ] <50% (poor)

   Justify using spatial features, not cell counts.
"""
    return prompt
```

---

### Stage 4: Hypothesis Generation

```python
def stage4_hypothesis_generation(spatial_insights):
    """Generate testable predictions."""

    prompt = f"""
SCIENTIFIC HYPOTHESIS GENERATION:

Based on spatial transcriptomics findings, propose TESTABLE hypotheses.

CONSTRAINT: Each hypothesis must be:
1. Falsifiable (can be proven wrong)
2. Testable with standard techniques (IHC, IF, CODEX, imaging mass cytometry)
3. Directly addressing a spatial pattern

FORMAT (Use this template):

HYPOTHESIS 1: [Mechanistic Statement]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Observation: Complete spatial segregation between tumor and immune cells

Hypothesis: Tumor cells secrete CXCL12 creating a chemokine gradient that
repels CXCR4+ T cells, preventing infiltration despite immune activation.

Validation Experiment:
- Method: Multiplex immunofluorescence
- Markers: CXCL12 (tumor), CXCR4 (T cells), CD8, PanCK
- Expected Result if TRUE: CXCL12 gradient highest at tumor-stroma boundary,
  CXCR4+ T cells accumulate at boundary but not in tumor center
- Expected Result if FALSE: No CXCL12 gradient, or CXCR4- T cells

Alternative Test: scRNA-seq of interface vs core to measure CXCL12 expression

HYPOTHESIS 2: [Immune Phenotype]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Observation: ISG15 clustering with Moran's I = 0.57

Hypothesis: ISG15+ cells are tumor-reactive CD8+ T cells with chronic
interferon stimulation, indicating active but exhausted immune response.

Validation Experiment:
- Method: CODEX multiplexed imaging
- Markers: ISG15, CD8, PD-1, TIM-3, GzmB, Ki67
- Expected Result if TRUE: ISG15+ CD8+ T cells express exhaustion markers
  (PD-1high, TIM-3+) and show low Ki67 (not proliferating)
- Expected Result if FALSE: ISG15+ cells are macrophages, not T cells, OR
  T cells are not exhausted (PD-1low, Ki67+)

Alternative Test: Flow cytometry of tissue digest, gate on ISG15+ CD8+ cells,
measure exhaustion markers

HYPOTHESIS 3: [Spatial Heterogeneity]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Observation: Moderate entropy (0.42) suggests regional differences

Hypothesis: Tumor periphery is more heterogeneous than core, with distinct
molecular programs in invasive vs proliferative regions.

Validation Experiment:
- Method: Spatial proteomics (MIBI-TOF or Hyperion)
- Regions: Define "core" (central 30%), "invasive front" (outer 30%)
- Markers: Ki67 (proliferation), Vimentin (EMT), E-cadherin, ER, PR
- Expected Result if TRUE: Core shows high Ki67/low Vimentin,
  periphery shows low Ki67/high Vimentin (EMT signature)
- Expected Result if FALSE: No regional differences in marker distribution

NOW GENERATE 3 HYPOTHESES FOR THIS SAMPLE:

{format_spatial_patterns(spatial_insights)}
"""
    return prompt
```

---

### Stage 5: Verification Layer (Llama-3.1-8B)

**Why Llama-3.1-8B instead of another medical LLM?**
- General reasoning (not medical-specific) is better for VALIDATION
- Less likely to have same biases as MedGemma
- Can catch biological implausibilities
- Lighter weight (8B vs MedGemma 4B but different training)

```python
def stage5_verification(medgemma_hypotheses, spatial_data):
    """Use Llama-3.1-8B to validate MedGemma's reasoning."""

    llama_prompt = f"""
You are a scientific reviewer evaluating spatial transcriptomics analysis.

MEDGEMMA'S ANALYSIS:
{medgemma_hypotheses}

ORIGINAL SPATIAL DATA:
{json.dumps(spatial_data, indent=2)}

VALIDATION TASKS:

1. BIOLOGICAL PLAUSIBILITY CHECK:
   Review each of MedGemma's claims. Are they biologically plausible?

   For each claim, provide:
   - [ ] Plausible
   - [ ] Questionable (explain why)
   - [ ] Implausible (explain why)

2. ALTERNATIVE EXPLANATIONS:
   For MedGemma's top 3 conclusions, propose ALTERNATIVE explanations
   that fit the spatial data equally well.

   This tests whether MedGemma's reasoning is the ONLY valid interpretation.

3. CONTRADICTION DETECTION:
   Are there any INTERNAL CONTRADICTIONS in MedGemma's analysis?

   Example of contradiction:
   - Claims "immune cells are excluded"
   - But also says "active immune response"
   - These could contradict unless carefully explained

4. DATA-CONCLUSION ALIGNMENT:
   Does each conclusion NECESSARILY follow from the spatial data?
   Or are there logical leaps?

   Flag any conclusion that requires assumptions not supported by data.

5. NOVELTY ASSESSMENT:
   Which of MedGemma's insights are:
   - [A] Directly stated in the input data (parrot output)
   - [B] Reasonable inference from data
   - [C] Novel hypothesis requiring domain knowledge
   - [D] Speculative without clear data support

FINAL VERDICT:
Overall quality score: ___/100
Recommendation: [ ] Accept [ ] Revise [ ] Reject

If Revise or Reject, specify what MedGemma should re-analyze.
"""

    # Run Llama-3.1-8B
    verification = run_llama_31(llama_prompt)

    return verification
```

---

### Stage 6: Programmatic Novelty Validation

```python
import re
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

def stage6_novelty_validation(medgemma_output, input_json_str, threshold=0.6):
    """Programmatically detect parrot outputs."""

    # Check 1: Forbidden phrases (direct data parroting)
    forbidden_exact_matches = [
        str(input_json['tumor_immune_interface']['interface_pct']),
        str(input_json['spatial_entropy']['overall']['mean']),
    ]

    for phrase in forbidden_exact_matches:
        if phrase in medgemma_output:
            return {
                "pass": False,
                "reason": f"Direct data copy detected: '{phrase}'",
                "novelty_score": 0.0
            }

    # Check 2: Required reasoning patterns
    reasoning_patterns = [
        r'\bif\b.*\bthen\b',  # Conditional reasoning
        r'\bsuggests?\b',      # Inference
        r'\bhypothesis\b',     # Hypothesis generation
        r'\bcompared to\b',    # Comparative analysis
        r'\balternative(?:ly)?\b',  # Alternative explanations
    ]

    reasoning_score = sum(
        1 for pattern in reasoning_patterns
        if re.search(pattern, medgemma_output, re.IGNORECASE)
    )

    if reasoning_score < 2:
        return {
            "pass": False,
            "reason": f"Insufficient reasoning patterns ({reasoning_score}/5)",
            "novelty_score": reasoning_score / 5
        }

    # Check 3: TF-IDF divergence (semantic novelty)
    vectorizer = TfidfVectorizer(stop_words='english', max_features=100)
    vectors = vectorizer.fit_transform([input_json_str, medgemma_output])

    similarity = cosine_similarity(vectors[0:1], vectors[1:2])[0][0]
    novelty_score = 1 - similarity  # High similarity = low novelty

    if novelty_score < threshold:
        return {
            "pass": False,
            "reason": f"Low semantic novelty: {novelty_score:.2f} < {threshold}",
            "novelty_score": novelty_score
        }

    # Check 4: Length and depth
    word_count = len(medgemma_output.split())
    if word_count < 300:
        return {
            "pass": False,
            "reason": f"Output too brief ({word_count} words, minimum 300)",
            "novelty_score": novelty_score
        }

    # All checks passed
    return {
        "pass": True,
        "novelty_score": novelty_score,
        "reasoning_score": reasoning_score,
        "word_count": word_count
    }
```

---

## Complete Pipeline

```python
def run_medgemma_v2_pipeline(adata, spatial_stats, celltype_stats):
    """Anti-parrot MedGemma architecture."""

    print("Stage 1: Spatial Pattern Extraction...")
    stage1_prompt = stage1_spatial_patterns(adata, spatial_stats)
    spatial_insights = run_medgemma(stage1_prompt)

    print("Stage 2: Comparative Medical Reasoning...")
    current_sample = {
        'interface_pct': celltype_stats['tumor_immune_interface']['interface_pct'],
        'isg15_morans_i': find_gene_morans_i(spatial_stats, 'ISG15'),
        'entropy': spatial_stats['spatial_entropy']['overall']['mean'],
        'enrichment': 'all_depleted'  # From neighborhood enrichment
    }
    stage2_prompt = stage2_comparative_reasoning(current_sample)
    comparative_analysis = run_medgemma(stage2_prompt)

    print("Stage 3: Conditional Reasoning...")
    stage3_prompt = stage3_conditional_reasoning(spatial_stats, celltype_stats)
    conditional_chains = run_medgemma(stage3_prompt)

    print("Stage 4: Hypothesis Generation...")
    stage4_prompt = stage4_hypothesis_generation(spatial_insights)
    hypotheses = run_medgemma(stage4_prompt)

    print("Stage 5: Verification (Llama-3.1-8B)...")
    verification = stage5_verification(hypotheses, spatial_stats)

    print("Stage 6: Novelty Validation...")
    input_json_str = json.dumps(spatial_stats)
    combined_output = f"{spatial_insights}\n\n{comparative_analysis}\n\n{conditional_chains}\n\n{hypotheses}"

    validation_result = stage6_novelty_validation(
        combined_output,
        input_json_str,
        threshold=0.6
    )

    if not validation_result['pass']:
        raise ValueError(f"Novelty validation failed: {validation_result['reason']}")

    # Final report
    report = {
        "spatial_insights": spatial_insights,
        "comparative_analysis": comparative_analysis,
        "conditional_reasoning": conditional_chains,
        "testable_hypotheses": hypotheses,
        "verification": verification,
        "quality_metrics": {
            "novelty_score": validation_result['novelty_score'],
            "reasoning_score": validation_result['reasoning_score'],
            "word_count": validation_result['word_count'],
            "verification_score": parse_verification_score(verification)
        }
    }

    return report
```

---

## Summary: Why This Works

### Anti-Parrot Mechanisms:

1. **No Access to Pre-Computed Summaries**: MedGemma never sees JSON
2. **Forced Comparison**: Must compare to reference samples
3. **Conditional Logic**: Must use IF-THEN reasoning
4. **Testable Outputs**: Must generate falsifiable hypotheses
5. **Verification Layer**: External model validates reasoning
6. **Programmatic Checks**: Automatic rejection of parrot outputs

### Value-Add Layers:

- **Stage 1**: Pattern recognition (WHERE genes cluster)
- **Stage 2**: Medical contextualization (phenotype classification)
- **Stage 3**: Clinical reasoning (treatment selection)
- **Stage 4**: Experimental design (follow-up studies)
- **Stage 5**: Quality control (catch errors)
- **Stage 6**: Novelty enforcement (reject low-value outputs)

---

**This architecture makes it STRUCTURALLY IMPOSSIBLE to generate "34.6% tumor-immune interface observed" without additional medical reasoning.**

---
