"""
MedGemma Self-Audit Stage
Detects fundamental errors, tangents, and fragility in generated reports.
"""

import re
import json
from typing import Dict, List, Tuple


# ==============================================================================
# SELF-AUDIT STAGE: FUNDAMENTAL ERROR DETECTION
# ==============================================================================

def detect_parroting(medgemma_output: str, input_features: Dict) -> Dict:
    """
    Detect if MedGemma is just parroting input data.

    Returns:
        audit_result: Dict with parroting detection results
    """
    violations = []

    # Extract numeric values from input
    def extract_numbers(obj, path=""):
        """Recursively extract all numeric values from nested dict."""
        numbers = []
        if isinstance(obj, dict):
            for k, v in obj.items():
                numbers.extend(extract_numbers(v, f"{path}.{k}"))
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                numbers.extend(extract_numbers(item, f"{path}[{i}]"))
        elif isinstance(obj, (int, float)):
            numbers.append((path, obj))
        return numbers

    input_numbers = extract_numbers(input_features)

    # Check if exact numeric values appear in output
    for path, value in input_numbers:
        # Format number as it might appear in text
        if isinstance(value, float):
            # Check multiple precision levels
            for precision in [0, 1, 2, 3]:
                formatted = f"{value:.{precision}f}"
                if formatted in medgemma_output:
                    violations.append({
                        'type': 'NUMERIC_PARROTING',
                        'path': path,
                        'value': value,
                        'formatted': formatted,
                        'severity': 'CRITICAL'
                    })
                    break
        elif str(value) in medgemma_output:
            violations.append({
                'type': 'NUMERIC_PARROTING',
                'path': path,
                'value': value,
                'severity': 'CRITICAL'
            })

    # Check for generic phrases (canonical laundering)
    generic_phrases = [
        "spatial segregation",
        "spatial organization",
        "spatial patterns",
        "observed for all cell types",
        "suggests",
        "indicates",
        "may indicate",
        "is characterized by"
    ]

    generic_count = 0
    for phrase in generic_phrases:
        count = medgemma_output.lower().count(phrase.lower())
        generic_count += count
        if count > 2:
            violations.append({
                'type': 'GENERIC_PHRASE_OVERUSE',
                'phrase': phrase,
                'count': count,
                'severity': 'WARNING'
            })

    if generic_count > 8:
        violations.append({
            'type': 'EXCESSIVE_GENERIC_LANGUAGE',
            'total_generic_phrases': generic_count,
            'severity': 'CRITICAL'
        })

    return {
        'parroting_detected': len([v for v in violations if v['severity'] == 'CRITICAL']) > 0,
        'violations': violations,
        'total_violations': len(violations)
    }


def detect_tangents(medgemma_output: str, expected_topics: List[str]) -> Dict:
    """
    Detect if MedGemma is going on tangents unrelated to spatial patterns.

    Returns:
        tangent_audit: Dict with tangent detection results
    """
    # Keywords that should appear in spatial transcriptomics report
    spatial_keywords = [
        'spatial', 'cluster', 'neighborhood', 'niche', 'interface',
        'heterogeneity', 'pattern', 'region', 'localization', 'proximity'
    ]

    # Keywords indicating potential tangents
    tangent_keywords = [
        'genomic instability', 'mutation', 'sequencing depth',
        'batch effect', 'normalization method', 'algorithm',
        'generally', 'typically', 'often', 'usually'
    ]

    output_lower = medgemma_output.lower()

    spatial_count = sum(output_lower.count(kw) for kw in spatial_keywords)
    tangent_count = sum(output_lower.count(kw) for kw in tangent_keywords)

    # Check sentence-level focus
    sentences = re.split(r'[.!?]+', medgemma_output)
    off_topic_sentences = []

    for i, sent in enumerate(sentences):
        sent_lower = sent.lower()
        sent_spatial = sum(kw in sent_lower for kw in spatial_keywords)
        sent_tangent = sum(kw in sent_lower for kw in tangent_keywords)

        if sent_tangent > 0 and sent_spatial == 0 and len(sent.split()) > 5:
            off_topic_sentences.append({
                'sentence_num': i,
                'text': sent.strip(),
                'tangent_keywords': [kw for kw in tangent_keywords if kw in sent_lower]
            })

    tangent_detected = tangent_count > spatial_count * 0.3  # More than 30% tangent vs spatial

    return {
        'tangent_detected': tangent_detected,
        'spatial_keyword_count': spatial_count,
        'tangent_keyword_count': tangent_count,
        'off_topic_sentences': off_topic_sentences,
        'tangent_ratio': tangent_count / max(spatial_count, 1)
    }


def detect_uncertainty_omission(medgemma_output: str, uncertainty_data: Dict) -> Dict:
    """
    Check if MedGemma mentions uncertainty, p-values, CIs when available.

    Returns:
        uncertainty_audit: Dict with uncertainty reporting assessment
    """
    # Keywords indicating uncertainty awareness
    uncertainty_keywords = [
        'p-value', 'p <', 'p=', 'confidence interval', 'ci',
        'uncertain', 'variability', 'variance', 'error',
        'low confidence', 'high confidence', 'moderate',
        'bootstrap', 'permutation', 'significant', 'non-significant'
    ]

    output_lower = medgemma_output.lower()
    uncertainty_mentions = sum(kw in output_lower for kw in uncertainty_keywords)

    # Check if input had uncertainty data
    has_pvalues = False
    has_ci = False

    if 'stage_1_spatial_patterns' in uncertainty_data:
        morans = uncertainty_data['stage_1_spatial_patterns'].get('morans_i', {})
        if 'p_value' in morans:
            has_pvalues = True
        if 'ci_lower' in morans:
            has_ci = True

    violations = []

    if has_pvalues and 'p-value' not in output_lower and 'p <' not in output_lower:
        violations.append({
            'type': 'MISSING_P_VALUE_REPORTING',
            'severity': 'CRITICAL',
            'message': 'P-values available in input but not reported'
        })

    if has_ci and 'confidence interval' not in output_lower and 'ci' not in output_lower:
        violations.append({
            'type': 'MISSING_CI_REPORTING',
            'severity': 'WARNING',
            'message': 'Confidence intervals available but not reported'
        })

    if uncertainty_mentions == 0 and (has_pvalues or has_ci):
        violations.append({
            'type': 'ZERO_UNCERTAINTY_LANGUAGE',
            'severity': 'CRITICAL',
            'message': 'No uncertainty language despite quantitative uncertainty data'
        })

    return {
        'uncertainty_omission_detected': len(violations) > 0,
        'uncertainty_keyword_count': uncertainty_mentions,
        'violations': violations
    }


def detect_overconfident_claims(medgemma_output: str) -> Dict:
    """
    Detect overconfident language that doesn't acknowledge limitations.

    Returns:
        overconfidence_audit: Dict with overconfidence detection
    """
    # Overconfident phrases
    overconfident_patterns = [
        r'\b(definitively|certainly|clearly|obviously|undoubtedly)\b',
        r'\b(proves?|demonstrates?|confirms?)\b',
        r'\b(will|must|cannot)\b',
        r'\b(always|never)\b'
    ]

    # Hedging phrases (good)
    hedging_patterns = [
        r'\b(may|might|could|possibly|potentially)\b',
        r'\b(suggests?|indicates?|appears?|seems?)\b',
        r'\b(likely|unlikely|probable)\b',
        r'\b(warrants? further|requires? validation)\b'
    ]

    overconfident_matches = []
    for pattern in overconfident_patterns:
        matches = re.finditer(pattern, medgemma_output, re.IGNORECASE)
        for match in matches:
            # Get context (20 chars before and after)
            start = max(0, match.start() - 20)
            end = min(len(medgemma_output), match.end() + 20)
            context = medgemma_output[start:end]
            overconfident_matches.append({
                'pattern': pattern,
                'match': match.group(),
                'context': context.strip()
            })

    hedging_count = sum(len(re.findall(p, medgemma_output, re.IGNORECASE)) for p in hedging_patterns)
    overconfident_count = len(overconfident_matches)

    overconfidence_detected = overconfident_count > hedging_count * 0.5

    return {
        'overconfidence_detected': overconfidence_detected,
        'overconfident_phrase_count': overconfident_count,
        'hedging_phrase_count': hedging_count,
        'overconfident_examples': overconfident_matches[:5]  # Top 5 examples
    }


def check_claim_evidence_alignment(medgemma_output: str, uncertainty_data: Dict) -> Dict:
    """
    Verify claims are supported by signal strength in data.

    Returns:
        alignment_audit: Dict with claim-evidence alignment assessment
    """
    violations = []

    # Extract signal strengths from uncertainty data
    signal_strengths = {}
    if 'stage_1_spatial_patterns' in uncertainty_data:
        morans = uncertainty_data['stage_1_spatial_patterns'].get('morans_i', {})
        if 'genes' in morans and 'signal_strength' in morans:
            for gene, strength in zip(morans['genes'], morans['signal_strength']):
                signal_strengths[gene] = strength

    # Check for weak signal genes being discussed as strong
    weak_genes = [g for g, s in signal_strengths.items() if s in ['WEAK', 'NONE']]
    strong_language_patterns = [
        r'(strongly|significantly|robustly) .{{0,30}}{gene}',
        r'{gene} .{{0,30}}(strongly|significantly|robustly)'
    ]

    for gene in weak_genes:
        for pattern_template in strong_language_patterns:
            pattern = pattern_template.replace('{gene}', re.escape(gene))
            if re.search(pattern, medgemma_output, re.IGNORECASE):
                violations.append({
                    'type': 'WEAK_SIGNAL_STRONG_CLAIM',
                    'gene': gene,
                    'actual_signal': signal_strengths[gene],
                    'severity': 'CRITICAL'
                })

    # Check if stopping decision is mentioned when signal is weak
    if 'stage_stopping' in uncertainty_data:
        decision = uncertainty_data['stage_stopping'].get('decision')
        if decision in ['STOP_WEAK_SIGNAL', 'STOP_LOW_QUALITY']:
            if 'weak signal' not in medgemma_output.lower() and 'low quality' not in medgemma_output.lower():
                violations.append({
                    'type': 'IGNORED_STOPPING_DECISION',
                    'decision': decision,
                    'severity': 'CRITICAL',
                    'message': 'Report should acknowledge weak signal and stop analysis'
                })

    return {
        'misalignment_detected': len(violations) > 0,
        'violations': violations
    }


# ==============================================================================
# FRAGILITY TESTING
# ==============================================================================

def test_prompt_fragility(medgemma_output_1: str, medgemma_output_2: str) -> Dict:
    """
    Compare two outputs from same data with slightly different prompts.
    High inconsistency = fragile reasoning.

    Returns:
        fragility_metrics: Dict with consistency assessment
    """
    # Tokenize into sentences
    sentences_1 = set(re.split(r'[.!?]+', medgemma_output_1.lower()))
    sentences_2 = set(re.split(r'[.!?]+', medgemma_output_2.lower()))

    # Jaccard similarity (sentence level)
    intersection = len(sentences_1 & sentences_2)
    union = len(sentences_1 | sentences_2)
    jaccard = intersection / max(union, 1)

    # Extract key claims (sentences with "suggests", "indicates", etc.)
    claim_keywords = ['suggests', 'indicates', 'shows', 'demonstrates', 'reveals']

    def extract_claims(text):
        claims = []
        for sent in re.split(r'[.!?]+', text):
            if any(kw in sent.lower() for kw in claim_keywords):
                claims.append(sent.strip().lower())
        return set(claims)

    claims_1 = extract_claims(medgemma_output_1)
    claims_2 = extract_claims(medgemma_output_2)

    claim_consistency = len(claims_1 & claims_2) / max(len(claims_1 | claims_2), 1)

    fragile = claim_consistency < 0.5  # <50% claim consistency = fragile

    return {
        'fragility_detected': fragile,
        'sentence_similarity': jaccard,
        'claim_consistency': claim_consistency,
        'unique_claims_output1': len(claims_1 - claims_2),
        'unique_claims_output2': len(claims_2 - claims_1)
    }


# ==============================================================================
# MAIN SELF-AUDIT FUNCTION
# ==============================================================================

def run_self_audit(medgemma_output: str, input_features: Dict, uncertainty_data: Dict) -> Dict:
    """
    Execute complete self-audit on MedGemma output.

    Returns:
        audit_report: Comprehensive audit results with PASS/FAIL decision
    """
    print("=" * 80)
    print("MEDGEMMA SELF-AUDIT STAGE")
    print("=" * 80)

    # Run all audit checks
    print("\n[1/6] Parroting detection...")
    parroting_audit = detect_parroting(medgemma_output, input_features)

    print("[2/6] Tangent detection...")
    tangent_audit = detect_tangents(medgemma_output, expected_topics=['spatial patterns'])

    print("[3/6] Uncertainty omission check...")
    uncertainty_audit = detect_uncertainty_omission(medgemma_output, uncertainty_data)

    print("[4/6] Overconfident claims detection...")
    overconfidence_audit = detect_overconfident_claims(medgemma_output)

    print("[5/6] Claim-evidence alignment...")
    alignment_audit = check_claim_evidence_alignment(medgemma_output, uncertainty_data)

    # Compile results
    audit_report = {
        'parroting': parroting_audit,
        'tangents': tangent_audit,
        'uncertainty_omission': uncertainty_audit,
        'overconfidence': overconfidence_audit,
        'claim_evidence_alignment': alignment_audit
    }

    # Overall decision
    critical_failures = []

    if parroting_audit['parroting_detected']:
        critical_failures.append("PARROTING_DETECTED")

    if uncertainty_audit['uncertainty_omission_detected']:
        critical_failures.append("UNCERTAINTY_OMISSION")

    if alignment_audit['misalignment_detected']:
        critical_failures.append("CLAIM_EVIDENCE_MISALIGNMENT")

    if tangent_audit['tangent_detected']:
        critical_failures.append("TANGENT_DETECTED")

    if overconfidence_audit['overconfidence_detected']:
        critical_failures.append("OVERCONFIDENCE_DETECTED")

    overall_pass = len(critical_failures) == 0

    audit_report['overall'] = {
        'pass': overall_pass,
        'critical_failures': critical_failures,
        'total_critical_failures': len(critical_failures)
    }

    print("\n" + "=" * 80)
    print("SELF-AUDIT RESULTS")
    print("=" * 80)
    print(f"Overall: {'PASS' if overall_pass else 'FAIL'}")
    if not overall_pass:
        print(f"Critical failures: {', '.join(critical_failures)}")
    print("=" * 80)

    return audit_report


if __name__ == "__main__":
    # Example usage
    example_output = """
    The tissue shows 34.6% tumor-immune interface. Spatial clustering is observed.
    ISG15, C1QA, C1QB show elevated expression. Neighborhood enrichment patterns
    are observed for all cell types, suggesting spatial segregation.
    """

    example_input = {
        'tumor_immune_interface': {'interface_pct': 34.6},
        'spatial_autocorrelation': {'top_genes': ['ISG15', 'C1QA', 'C1QB']}
    }

    example_uncertainty = {
        'stage_1_spatial_patterns': {
            'morans_i': {
                'genes': ['ISG15', 'C1QA'],
                'p_value': [0.001, 0.045],
                'signal_strength': ['STRONG', 'WEAK']
            }
        }
    }

    audit = run_self_audit(example_output, example_input, example_uncertainty)

    print("\nExample Audit Results:")
    print(json.dumps(audit, indent=2))
