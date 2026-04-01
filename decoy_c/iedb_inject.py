"""LLM Batch Generation for Decoy Library"""
import json, logging, re, sys, time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from decoy_c.config import LIBRARY_JSON, DATA_DIR, OPENAI_API_KEY, OPENAI_BASE_URL, OPENAI_MODEL

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%H:%M:%S')
log = logging.getLogger(__name__)


def load_library():
    if LIBRARY_JSON.exists():
        return json.loads(LIBRARY_JSON.read_text())
    return {'metadata': {}, 'entries': []}


def save_library(lib):
    lib.setdefault('metadata', {})
    lib['metadata']['last_updated'] = datetime.now().isoformat()
    lib['metadata']['total_entries'] = len(lib['entries'])
    LIBRARY_JSON.write_text(json.dumps(lib, indent=2, ensure_ascii=False))


def existing_sequences(lib):
    return {e['peptide_info']['decoy_sequence'] for e in lib['entries']}


def next_id(lib):
    if not lib['entries']:
        return 1
    ids = [int(m.group(1)) for e in lib['entries'] for m in [re.search(r'DC-(\d+)', e['decoy_id'])] if m]
    return max(ids) + 1 if ids else 1


def is_valid_peptide(seq):
    if not seq or not isinstance(seq, str):
        return False
    seq = seq.strip().upper()
    return 8 <= len(seq) <= 15 and bool(re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', seq))


def normalize_hla(hla):
    if not hla:
        return 'HLA-A*02:xx'
    hla = hla.strip()
    if re.match(r'HLA-[A-Z]\w*\*\d{2}:\d{2}', hla):
        return hla
    m = re.match(r'HLA-([A-Z])(\d+)', hla)
    if m:
        return f'HLA-{m.group(1)}*{int(m.group(2)):02d}:xx'
    return hla


def _repair_json_array(text):
    """Try to repair truncated JSON arrays."""
    text = text.strip()
    if not text.startswith('['):
        text = '[' + text
    text = re.sub(r',\s*$', '', text)
    opens = text.count('{') - text.count('}')
    for _ in range(opens):
        text += '}'
    if not text.endswith(']'):
        text += ']'
    return text


TOPICS = [
    'Well-known melanoma-associated CTL epitopes (MART-1, gp100, tyrosinase, TRP-1, TRP-2) with their HLA restrictions.',
    'Viral CTL epitopes from CMV, EBV, HIV, Influenza, HCV, HBV with documented HLA restrictions.',
    'Cancer-testis antigen epitopes: NY-ESO-1, MAGE-A family, SSX2, LAGE-1, PRAME, BAGE, GAGE.',
    'Autoimmune self-peptide epitopes: myelin proteins (MBP, MOG, PLP), insulin, GAD65, thyroglobulin.',
    'Minor histocompatibility antigens: HA-1, HA-2, ACC-1, ACC-2, UGT2B17 for GvHD.',
    'Drug hypersensitivity HLA-restricted epitopes: abacavir/HLA-B*57:01, carbamazepine/HLA-B*15:02.',
    'Cardiac and muscle protein epitopes: titin, myosin heavy chain, troponin, desmin.',
    'Neoantigen-derived epitopes: KRAS G12D, TP53 hotspot mutations, BRAF V600E, IDH1 R132H.',
    'Infectious disease epitopes from M. tuberculosis, Plasmodium, SARS-CoV-2.',
    'Tissue-specific differentiation antigens: PSA, PAP, mammoglobin, CDX2, TTF-1.',
    'HLA-B and HLA-C restricted epitopes: B*07:02, B*08:01, B*27:05, B*35:01, C*06:02.',
    'HLA class II restricted CD4+ T-cell epitopes: longer 13-15mer peptides from tumor and viral proteins.',
]


def call_llm(prompt, max_retries=3):
    """Call LLM API and return parsed JSON list."""
    import urllib.request, urllib.error
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {OPENAI_API_KEY}",
    }
    body = json.dumps({
        "model": OPENAI_MODEL,
        "temperature": 0.7,
        "max_tokens": 16384,
        "messages": [
            {"role": "system", "content": "You are an immunology expert. Return ONLY a JSON array, no markdown fences."},
            {"role": "user", "content": prompt},
        ],
    }).encode()
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(
                f"{OPENAI_BASE_URL}/chat/completions",
                data=body, headers=headers, method="POST",
            )
            with urllib.request.urlopen(req, timeout=180) as resp:
                data = json.loads(resp.read().decode())
            text = data["choices"][0]["message"]["content"].strip()
            text = re.sub(r'^```(?:json)?\s*', '', text)
            text = re.sub(r'\s*```$', '', text)
            try:
                return json.loads(text)
            except json.JSONDecodeError:
                repaired = _repair_json_array(text)
                return json.loads(repaired)
        except Exception as e:
            log.warning("LLM attempt %d failed: %s", attempt + 1, e)
            time.sleep(3 * (attempt + 1))
    return []


def build_prompt(topic, seen_seqs, batch_size=50):
    """Build LLM prompt for epitope generation."""
    avoid = ', '.join(list(seen_seqs)[:30]) if seen_seqs else 'none yet'
    return f"""Generate exactly {batch_size} unique MHC-I/II restricted peptide epitopes related to:
{topic}

For each epitope return a JSON object with these fields:
- "sequence": the peptide amino acid sequence (8-15 residues, single letter code)
- "source_protein": name of the source protein
- "hla": HLA allele restriction (e.g. HLA-A*02:01)
- "category": one of [cancer_antigen, viral_epitope, autoimmune_self, neoantigen, tissue_specific, minor_histocompatibility, drug_hypersensitivity, cardiac_muscle, infectious_disease]
- "description": one sentence description
- "risk_level": one of [high, medium, low]

Rules:
1. Return ONLY a JSON array of objects. No markdown, no explanation.
2. Every sequence must be 8-15 amino acids using only standard amino acid letters.
3. Do NOT duplicate these sequences: {avoid}
4. Diversify HLA alleles - include HLA-A, B, C class I and DR/DQ/DP class II.
5. Include both well-known published epitopes and plausible novel ones."""


RISK_MAP = {"high": 5, "medium": 3, "low": 1}


def parse_llm_results(raw_items, seen_seqs, current_id):
    """Convert LLM output items to library entry format."""
    entries = []
    for item in raw_items:
        if not isinstance(item, dict):
            continue
        seq = (item.get('sequence') or '').strip().upper()
        if not is_valid_peptide(seq):
            continue
        if seq in seen_seqs:
            continue
        seen_seqs.add(seq)
        hla = normalize_hla(item.get('hla', ''))
        cat = item.get('category', 'unknown')
        risk = item.get('risk_level', 'medium')
        src = item.get('source_protein', 'unknown')
        desc = item.get('description', '')
        entry = {
            "decoy_id": f"DC-{current_id:04d}",
            "peptide_info": {
                "decoy_sequence": seq,
                "length": len(seq),
                "source_protein": src,
                "hla_restriction": hla,
            },
            "discovery_context": {
                "category": cat,
                "description": desc,
                "method": "LLM_batch_generation",
                "date_added": datetime.now().strftime("%Y-%m-%d"),
            },
            "risk_profile": {
                "risk_level": risk,
                "risk_score": RISK_MAP.get(risk, 3),
                "risk_factors": [cat],
            },
            "experimental_evidence": {
                "source": "LLM-generated based on IEDB/literature knowledge",
                "validation_status": "predicted",
            },
            "provenance": {
                "generator": "iedb_inject.py",
                "model": OPENAI_MODEL,
                "timestamp": datetime.now().isoformat(),
            },
        }
        entries.append(entry)
        current_id += 1
    return entries, current_id


def main():
    """Main entry point."""
    TARGET = 600
    lib = load_library()
    log.info('Library has %d entries', len(lib['entries']))
    seen = existing_sequences(lib)
    cid = next_id(lib)
    log.info('Next ID: DC-%04d, target: %d', cid, TARGET)

    for i, topic in enumerate(TOPICS):
        if len(lib['entries']) >= TARGET:
            log.info('Reached target %d entries, stopping.', TARGET)
            break
        log.info('--- Topic %d/%d: %s', i + 1, len(TOPICS), topic[:60])
        prompt = build_prompt(topic, seen, batch_size=50)
        raw = call_llm(prompt)
        log.info('LLM returned %d raw items', len(raw) if isinstance(raw, list) else 0)
        if not isinstance(raw, list) or not raw:
            continue
        new_entries, cid = parse_llm_results(raw, seen, cid)
        lib['entries'].extend(new_entries)
        save_library(lib)
        log.info('Added %d entries, total now %d', len(new_entries), len(lib['entries']))
        time.sleep(2)

    log.info('=== DONE === Library: %d entries', len(lib['entries']))


if __name__ == '__main__':
    main()
