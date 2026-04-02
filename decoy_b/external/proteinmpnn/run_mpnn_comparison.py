#!/usr/bin/env python
"""
Run ProteinMPNN on Tfold-generated PDB structures to predict sequences
at CDR3 (masked) positions only, and save results in matching format
for comparison with Tfold's predictions.

Usage:
    conda activate S3AI_new
    python run_mpnn_comparison.py

Output:
    - output_results/<sample_name>_probs.npz
    - output_results/<sample_name>_probs.json
    - comparison_report.md
"""

import os
import sys
import json
import glob
import copy
import time
import numpy as np
import torch
import torch.nn.functional as F

MPNN_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, MPNN_DIR)

from protein_mpnn_utils import (
    parse_PDB, tied_featurize, StructureDatasetPDB, ProteinMPNN, _S_to_seq
)

# Configuration
TFOLD_OUTPUT_DIR = "/share/liuyutian/tfold/batch_cdr3_generation/output_results"
MPNN_OUTPUT_DIR = os.path.join(MPNN_DIR, "output_results")
MODEL_WEIGHTS_DIR = os.path.join(MPNN_DIR, "vanilla_model_weights")
MODEL_NAME = "v_48_020"
SEED = 42
ALPHABET_20 = "ACDEFGHIKLMNPQRSTVWY"
ALPHABET_21 = "ACDEFGHIKLMNPQRSTVWYX"


def load_tfold_data(sample_name):
    json_path = os.path.join(TFOLD_OUTPUT_DIR, f"{sample_name}_probs.json")
    npz_path = os.path.join(TFOLD_OUTPUT_DIR, f"{sample_name}_probs.npz")
    with open(json_path, 'r') as f:
        json_data = json.load(f)
    npz_data = np.load(npz_path, allow_pickle=True)
    return json_data, npz_data


def get_mask_positions_from_tfold(sample_name):
    npz_path = os.path.join(TFOLD_OUTPUT_DIR, f"{sample_name}_probs.npz")
    npz_data = np.load(npz_path, allow_pickle=True)
    mask = npz_data['mask']
    masked_positions = np.where(mask == 1)[0].tolist()
    sequence = str(npz_data['sequence'])
    return masked_positions, mask, sequence


def build_fixed_positions(pdb_dict, masked_positions_0indexed, chain_order):
    chain_lengths = []
    chain_names = []
    for ch in chain_order:
        key = f"seq_chain_{ch}"
        if key in pdb_dict:
            chain_lengths.append(len(pdb_dict[key]))
            chain_names.append(ch)

    designed_by_chain = {ch: [] for ch in chain_names}
    offset = 0
    for i, ch in enumerate(chain_names):
        ch_len = chain_lengths[i]
        for pos in masked_positions_0indexed:
            if offset <= pos < offset + ch_len:
                pos_in_chain_1indexed = pos - offset + 1
                designed_by_chain[ch].append(pos_in_chain_1indexed)
        offset += ch_len

    fixed_positions = {}
    offset = 0
    for i, ch in enumerate(chain_names):
        ch_len = chain_lengths[i]
        all_positions = set(range(1, ch_len + 1))
        designed_set = set(designed_by_chain[ch])
        fixed_set = sorted(all_positions - designed_set)
        fixed_positions[ch] = fixed_set
        offset += ch_len

    return fixed_positions


def run_mpnn_on_sample(model, sample_name, device):
    pdb_path = os.path.join(TFOLD_OUTPUT_DIR, f"{sample_name}.pdb")
    masked_positions, tfold_mask, tfold_sequence = get_mask_positions_from_tfold(sample_name)

    print(f"\nProcessing {sample_name}:")
    print(f"  Sequence length: {len(tfold_sequence)}")
    print(f"  Masked (CDR3) positions: {masked_positions}")

    pdb_dict_list = parse_PDB(pdb_path, ca_only=False)
    pdb_dict = pdb_dict_list[0]
    all_chain_list = [item[-1:] for item in list(pdb_dict) if item[:9] == 'seq_chain']
    print(f"  Chains in PDB: {all_chain_list}")

    # Build concatenated sequence from PDB dict (alphabetical chain order)
    concat_seq = ""
    chain_order_in_concat = []
    for ch in all_chain_list:
        key = f"seq_chain_{ch}"
        if key in pdb_dict:
            chain_order_in_concat.append(ch)
            concat_seq += pdb_dict[key]

    print(f"  MPNN concat sequence length: {len(concat_seq)}")

    # Get Tfold chain order from PDB REMARKs
    tfold_json, _ = load_tfold_data(sample_name)
    tfold_chain_order = []
    tfold_chain_seqs = {}
    with open(pdb_path, 'r') as f:
        for line in f:
            if "Predicted Sequence for chain" in line:
                parts = line.strip().split("chain ")
                if len(parts) >= 2:
                    rest = parts[1]
                    chain_id = rest[0]
                    seq = rest.split(": ")[1].strip()
                    tfold_chain_order.append(chain_id)
                    tfold_chain_seqs[chain_id] = seq
    print(f"  Tfold chain order: {tfold_chain_order}")

    # Build position mapping: tfold_pos -> (chain, pos_in_chain)
    tfold_offset = 0
    tfold_pos_to_chain = {}
    for ch in tfold_chain_order:
        ch_len = len(tfold_chain_seqs[ch])
        for i in range(ch_len):
            tfold_pos_to_chain[tfold_offset + i] = (ch, i)
        tfold_offset += ch_len

    # Build MPNN position mapping
    mpnn_offset = 0
    chain_to_mpnn_offset = {}
    for ch in chain_order_in_concat:
        key = f"seq_chain_{ch}"
        chain_to_mpnn_offset[ch] = mpnn_offset
        mpnn_offset += len(pdb_dict[key])

    # Map masked positions from Tfold to MPNN
    mpnn_masked_positions = []
    for tfold_pos in masked_positions:
        if tfold_pos in tfold_pos_to_chain:
            ch, pos_in_chain = tfold_pos_to_chain[tfold_pos]
            if ch in chain_to_mpnn_offset:
                mpnn_pos = chain_to_mpnn_offset[ch] + pos_in_chain
                mpnn_masked_positions.append(mpnn_pos)
    print(f"  MPNN masked positions: {mpnn_masked_positions}")

    # Build fixed positions
    fixed_positions = build_fixed_positions(pdb_dict, mpnn_masked_positions, chain_order_in_concat)
    chain_id_dict = {pdb_dict['name']: (all_chain_list, [])}
    fixed_positions_dict = {pdb_dict['name']: fixed_positions}

    designed_chain_list = list(set([tfold_pos_to_chain[p][0] for p in masked_positions if p in tfold_pos_to_chain]))
    print(f"  Designed chains: {designed_chain_list}")
    for ch in all_chain_list:
        if ch in fixed_positions:
            ch_len = len(pdb_dict.get(f"seq_chain_{ch}", ""))
            print(f"    Chain {ch}: {len(fixed_positions[ch])}/{ch_len} fixed")

    # Featurize
    batch_clones = [copy.deepcopy(pdb_dict_list[0])]
    X, S, mask, lengths, chain_M, chain_encoding_all, chain_list_list, \
        visible_list_list, masked_list_list, masked_chain_length_list_list, \
        chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask, \
        tied_pos_list_of_lists_list, pssm_coef, pssm_bias, \
        pssm_log_odds_all, bias_by_res_all, tied_beta = tied_featurize(
            batch_clones, device, chain_id_dict, fixed_positions_dict,
            omit_AA_dict=None, tied_positions_dict=None, pssm_dict=None,
            bias_by_res_dict=None, ca_only=False
        )

    total_length = X.shape[1]
    print(f"  Total MPNN length: {total_length}")

    # Conditional probs: p(s_i | rest of sequence and backbone)
    randn_1 = torch.randn(chain_M.shape, device=device)
    log_conditional_probs = model.conditional_probs(
        X, S, mask, chain_M * chain_M_pos, residue_idx, chain_encoding_all, randn_1,
        backbone_only=False
    )
    conditional_probs = torch.exp(log_conditional_probs)

    # Sampling
    omit_AAs_np = np.array([AA in 'X' for AA in ALPHABET_21]).astype(np.float32)
    bias_AAs_np = np.zeros(len(ALPHABET_21))
    randn_2 = torch.randn(chain_M.shape, device=device)
    sample_dict = model.sample(
        X, randn_2, S, chain_M, chain_encoding_all, residue_idx,
        mask=mask, temperature=0.1,
        omit_AAs_np=omit_AAs_np, bias_AAs_np=bias_AAs_np,
        chain_M_pos=chain_M_pos, omit_AA_mask=omit_AA_mask,
        pssm_coef=pssm_coef, pssm_bias=pssm_bias,
        pssm_multi=0.0, pssm_log_odds_flag=False,
        pssm_log_odds_mask=None, pssm_bias_flag=False,
        bias_by_res=bias_by_res_all
    )
    S_sample = sample_dict["S"]
    sample_probs = sample_dict["probs"]

    cond_probs_np = conditional_probs[0].cpu().numpy()  # [L, 21]
    sample_probs_np = sample_probs[0].cpu().numpy()  # [L, 21]

    # Build output probs (20 AAs, renormalized)
    output_probs = np.zeros((total_length, 20), dtype=np.float32)
    for i in range(min(total_length, cond_probs_np.shape[0])):
        p20 = cond_probs_np[i, :20]
        p_sum = p20.sum()
        if p_sum > 0:
            output_probs[i] = p20 / p_sum

    # Build mask in Tfold sequence space
    output_mask = np.zeros(len(tfold_sequence), dtype=np.int8)
    for pos in masked_positions:
        output_mask[pos] = 1

    # Build designed sequence
    mpnn_sequence_list = list(tfold_sequence)
    for tfold_pos in masked_positions:
        if tfold_pos in tfold_pos_to_chain:
            ch, pos_in_chain = tfold_pos_to_chain[tfold_pos]
            if ch in chain_to_mpnn_offset:
                mpnn_pos = chain_to_mpnn_offset[ch] + pos_in_chain
                if mpnn_pos < S_sample.shape[1]:
                    aa_idx = S_sample[0, mpnn_pos].item()
                    mpnn_sequence_list[tfold_pos] = ALPHABET_21[aa_idx]
    mpnn_sequence = "".join(mpnn_sequence_list)

    # Map probs to Tfold sequence space
    output_probs_tfold = np.zeros((len(tfold_sequence), 20), dtype=np.float32)
    for tfold_pos in range(len(tfold_sequence)):
        if tfold_pos in tfold_pos_to_chain:
            ch, pos_in_chain = tfold_pos_to_chain[tfold_pos]
            if ch in chain_to_mpnn_offset:
                mpnn_pos = chain_to_mpnn_offset[ch] + pos_in_chain
                if mpnn_pos < output_probs.shape[0]:
                    output_probs_tfold[tfold_pos] = output_probs[mpnn_pos]

    # Full PDB probs
    total_pdb_residues = sum(len(pdb_dict.get(f"seq_chain_{ch}", "")) for ch in all_chain_list)
    output_probs_full = np.zeros((total_pdb_residues, 20), dtype=np.float32)
    pdb_offset = 0
    for ch in chain_order_in_concat:
        key = f"seq_chain_{ch}"
        ch_len = len(pdb_dict.get(key, ""))
        mpnn_ch_offset = chain_to_mpnn_offset[ch]
        for i in range(ch_len):
            mpnn_pos = mpnn_ch_offset + i
            if mpnn_pos < output_probs.shape[0]:
                output_probs_full[pdb_offset + i] = output_probs[mpnn_pos]
        pdb_offset += ch_len

    # Build JSON output
    json_output = {
        "id": sample_name,
        "full_sequence": mpnn_sequence,
    }
    for ch in tfold_chain_order:
        chain_key = f"{ch}_chain_predicted_probs"
        chain_preds = []
        for tfold_pos in masked_positions:
            if tfold_pos in tfold_pos_to_chain:
                pos_ch, pos_in_chain = tfold_pos_to_chain[tfold_pos]
                if pos_ch == ch:
                    probs_20 = output_probs_tfold[tfold_pos]
                    predicted_aa_idx = np.argmax(probs_20)
                    predicted_aa = ALPHABET_20[predicted_aa_idx]
                    probs_dict = {aa: float(probs_20[j]) for j, aa in enumerate(ALPHABET_20)}
                    chain_preds.append({
                        "pos": pos_in_chain + 1,
                        "predicted_aa": predicted_aa,
                        "probs": probs_dict
                    })
        json_output[chain_key] = chain_preds

    amino_acids = np.array(list(ALPHABET_20))

    return {
        "json_output": json_output,
        "probs_full": output_probs_full,
        "probs_tfold_seq": output_probs_tfold,
        "mask": output_mask,
        "sequence": mpnn_sequence,
        "amino_acids": amino_acids,
        "masked_positions": masked_positions,
        "mpnn_masked_positions": mpnn_masked_positions,
    }


def compute_comparison_metrics(sample_name, mpnn_result):
    tfold_json, tfold_npz = load_tfold_data(sample_name)
    tfold_sequence = str(tfold_npz['sequence'])
    mpnn_sequence = mpnn_result['sequence']
    masked_positions = mpnn_result['masked_positions']

    metrics = {
        "sample_name": sample_name,
        "sequence_length": len(tfold_sequence),
        "num_masked_positions": len(masked_positions),
        "masked_positions": masked_positions,
    }

    tfold_masked_aas = [tfold_sequence[p] for p in masked_positions]
    mpnn_masked_aas = [mpnn_sequence[p] for p in masked_positions]
    metrics["tfold_masked_aas"] = "".join(tfold_masked_aas)
    metrics["mpnn_masked_aas"] = "".join(mpnn_masked_aas)

    matches = sum(1 for a, b in zip(tfold_masked_aas, mpnn_masked_aas) if a == b)
    metrics["agreement_rate"] = matches / len(masked_positions) if masked_positions else 0.0

    mpnn_json = mpnn_result['json_output']
    kl_divs = []
    cosine_sims = []
    tfold_entropies = []
    mpnn_entropies = []
    tfold_max_probs = []
    mpnn_max_probs = []
    chain_metrics = {}

    for chain_key in [k for k in tfold_json.keys() if k.endswith('_chain_predicted_probs')]:
        chain_id = chain_key.split('_')[0]
        tfold_chain_preds = tfold_json[chain_key]
        mpnn_chain_preds = mpnn_json.get(chain_key, [])
        if not mpnn_chain_preds:
            continue

        chain_kl_divs = []
        chain_cosine_sims = []
        chain_comparison = []

        for t_pred, m_pred in zip(tfold_chain_preds, mpnn_chain_preds):
            t_probs = np.array([t_pred['probs'].get(aa, 0.0) for aa in ALPHABET_20])
            m_probs = np.array([m_pred['probs'].get(aa, 0.0) for aa in ALPHABET_20])
            t_probs = t_probs / (t_probs.sum() + 1e-10)
            m_probs = m_probs / (m_probs.sum() + 1e-10)

            kl = np.sum(t_probs * np.log((t_probs + 1e-10) / (m_probs + 1e-10)))
            chain_kl_divs.append(kl)

            cos_sim = np.dot(t_probs, m_probs) / (np.linalg.norm(t_probs) * np.linalg.norm(m_probs) + 1e-10)
            chain_cosine_sims.append(cos_sim)

            t_entropy = -np.sum(t_probs * np.log(t_probs + 1e-10))
            m_entropy = -np.sum(m_probs * np.log(m_probs + 1e-10))
            tfold_entropies.append(t_entropy)
            mpnn_entropies.append(m_entropy)
            tfold_max_probs.append(float(t_probs.max()))
            mpnn_max_probs.append(float(m_probs.max()))

            chain_comparison.append({
                "pos": t_pred['pos'],
                "tfold_aa": t_pred['predicted_aa'],
                "mpnn_aa": m_pred['predicted_aa'],
                "match": t_pred['predicted_aa'] == m_pred['predicted_aa'],
                "kl_divergence": float(kl),
                "cosine_similarity": float(cos_sim),
                "tfold_confidence": float(t_probs.max()),
                "mpnn_confidence": float(m_probs.max()),
                "tfold_entropy": float(t_entropy),
                "mpnn_entropy": float(m_entropy),
            })

        chain_metrics[chain_id] = {
            "comparisons": chain_comparison,
            "mean_kl_divergence": float(np.mean(chain_kl_divs)) if chain_kl_divs else 0,
            "mean_cosine_similarity": float(np.mean(chain_cosine_sims)) if chain_cosine_sims else 0,
        }
        kl_divs.extend(chain_kl_divs)
        cosine_sims.extend(chain_cosine_sims)

    metrics["chain_metrics"] = chain_metrics
    metrics["overall_mean_kl_divergence"] = float(np.mean(kl_divs)) if kl_divs else 0
    metrics["overall_mean_cosine_similarity"] = float(np.mean(cosine_sims)) if cosine_sims else 0
    metrics["tfold_mean_entropy"] = float(np.mean(tfold_entropies)) if tfold_entropies else 0
    metrics["mpnn_mean_entropy"] = float(np.mean(mpnn_entropies)) if mpnn_entropies else 0
    metrics["tfold_mean_max_prob"] = float(np.mean(tfold_max_probs)) if tfold_max_probs else 0
    metrics["mpnn_mean_max_prob"] = float(np.mean(mpnn_max_probs)) if mpnn_max_probs else 0

    return metrics


def generate_report(all_metrics):
    report = []
    report.append("# Tfold vs ProteinMPNN: CDR3 Sequence Design Comparison Report\n")
    report.append(f"**Generated:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    report.append(f"**Number of samples:** {len(all_metrics)}\n")
    report.append(f"**ProteinMPNN model:** {MODEL_NAME} (vanilla)\n")
    report.append(f"**Temperature:** 0.1 (for sampling)\n")
    report.append(f"**Seed:** {SEED}\n")

    report.append("\n## Summary Statistics\n")

    agreement_rates = [m['agreement_rate'] for m in all_metrics]
    kl_divs = [m['overall_mean_kl_divergence'] for m in all_metrics]
    cos_sims = [m['overall_mean_cosine_similarity'] for m in all_metrics]
    tfold_entropies = [m['tfold_mean_entropy'] for m in all_metrics]
    mpnn_entropies = [m['mpnn_mean_entropy'] for m in all_metrics]
    tfold_max_probs = [m['tfold_mean_max_prob'] for m in all_metrics]
    mpnn_max_probs = [m['mpnn_mean_max_prob'] for m in all_metrics]

    report.append("| Metric | Mean | Std | Min | Max |")
    report.append("|--------|------|-----|-----|-----|")
    for name, vals in [
        ("Agreement Rate", agreement_rates),
        ("KL Divergence (Tfold||MPNN)", kl_divs),
        ("Cosine Similarity", cos_sims),
        ("Tfold Mean Entropy", tfold_entropies),
        ("MPNN Mean Entropy", mpnn_entropies),
        ("Tfold Mean Max Prob", tfold_max_probs),
        ("MPNN Mean Max Prob", mpnn_max_probs),
    ]:
        arr = np.array(vals)
        report.append(f"| {name} | {arr.mean():.4f} | {arr.std():.4f} | {arr.min():.4f} | {arr.max():.4f} |")

    report.append("\n\n## Interpretation\n")
    report.append("- **Agreement Rate**: Fraction of CDR3 positions where Tfold and ProteinMPNN predict the same amino acid.")
    report.append("- **KL Divergence**: Measures how different the probability distributions are (lower = more similar).")
    report.append("- **Cosine Similarity**: Measures similarity of probability vectors (higher = more similar).")
    report.append("- **Entropy**: Measures uncertainty/diversity of predictions (higher = more uncertain).")
    report.append("- **Max Prob**: Confidence of the top prediction (higher = more confident).\n")

    # Per-sample details
    report.append("\n## Per-Sample Results\n")
    report.append("| Sample | Seq Len | #Masked | Agreement | KL Div | Cos Sim | Tfold AAs | MPNN AAs |")
    report.append("|--------|---------|---------|-----------|--------|---------|-----------|----------|")
    for m in all_metrics:
        report.append(
            f"| {m['sample_name']} | {m['sequence_length']} | {m['num_masked_positions']} | "
            f"{m['agreement_rate']:.3f} | {m['overall_mean_kl_divergence']:.3f} | "
            f"{m['overall_mean_cosine_similarity']:.3f} | {m['tfold_masked_aas']} | {m['mpnn_masked_aas']} |"
        )

    # Per-position detail for first few samples
    report.append("\n\n## Detailed Per-Position Comparison (First 5 Samples)\n")
    for m in all_metrics[:5]:
        report.append(f"\n### {m['sample_name']}\n")
        report.append(f"- Tfold CDR3: `{m['tfold_masked_aas']}`")
        report.append(f"- MPNN  CDR3: `{m['mpnn_masked_aas']}`")
        report.append(f"- Agreement: {m['agreement_rate']:.3f}\n")

        for chain_id, chain_data in m['chain_metrics'].items():
            report.append(f"\n**Chain {chain_id}:**\n")
            report.append("| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |")
            report.append("|-----|----------|---------|-------|--------|---------|------------|-----------|")
            for c in chain_data['comparisons']:
                match_str = "✓" if c['match'] else "✗"
                report.append(
                    f"| {c['pos']} | {c['tfold_aa']} | {c['mpnn_aa']} | {match_str} | "
                    f"{c['kl_divergence']:.3f} | {c['cosine_similarity']:.3f} | "
                    f"{c['tfold_confidence']:.3f} | {c['mpnn_confidence']:.3f} |"
                )

    return "\n".join(report)


def main():
    print("=" * 60)
    print("Tfold vs ProteinMPNN CDR3 Comparison")
    print("=" * 60)

    # Set seeds
    torch.manual_seed(SEED)
    np.random.seed(SEED)

    # Create output directory
    os.makedirs(MPNN_OUTPUT_DIR, exist_ok=True)

    # Load model
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    checkpoint_path = os.path.join(MODEL_WEIGHTS_DIR, f"{MODEL_NAME}.pt")
    checkpoint = torch.load(checkpoint_path, map_location=device)
    print(f"Loaded model: {MODEL_NAME}, noise level: {checkpoint['noise_level']}A, edges: {checkpoint['num_edges']}")

    model = ProteinMPNN(
        ca_only=False, num_letters=21, node_features=128, edge_features=128,
        hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
        augment_eps=0.0, k_neighbors=checkpoint['num_edges']
    )
    model.to(device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    # Find all sample PDB files
    pdb_files = sorted(glob.glob(os.path.join(TFOLD_OUTPUT_DIR, "*.pdb")))
    sample_names = [os.path.basename(f).replace(".pdb", "") for f in pdb_files]
    print(f"\nFound {len(sample_names)} samples")

    all_metrics = []
    failed_samples = []

    with torch.no_grad():
        for idx, sample_name in enumerate(sample_names):
            try:
                print(f"\n{'='*40} [{idx+1}/{len(sample_names)}] {'='*40}")
                result = run_mpnn_on_sample(model, sample_name, device)

                # Save NPZ
                npz_path = os.path.join(MPNN_OUTPUT_DIR, f"{sample_name}_probs.npz")
                np.savez(
                    npz_path,
                    probs=result['probs_full'],
                    mask=result['mask'],
                    sequence=result['sequence'],
                    amino_acids=result['amino_acids'],
                )
                print(f"  Saved: {npz_path}")

                # Save JSON
                json_path = os.path.join(MPNN_OUTPUT_DIR, f"{sample_name}_probs.json")
                with open(json_path, 'w') as f:
                    json.dump(result['json_output'], f, indent=2)
                print(f"  Saved: {json_path}")

                # Compute comparison metrics
                metrics = compute_comparison_metrics(sample_name, result)
                all_metrics.append(metrics)

                print(f"  Agreement rate: {metrics['agreement_rate']:.3f}")
                print(f"  Tfold CDR3: {metrics['tfold_masked_aas']}")
                print(f"  MPNN  CDR3: {metrics['mpnn_masked_aas']}")
                print(f"  KL Div: {metrics['overall_mean_kl_divergence']:.3f}")
                print(f"  Cos Sim: {metrics['overall_mean_cosine_similarity']:.3f}")

            except Exception as e:
                print(f"  ERROR: {e}")
                import traceback
                traceback.print_exc()
                failed_samples.append((sample_name, str(e)))

    # Generate report
    print(f"\n{'='*60}")
    print(f"Generating comparison report...")
    report = generate_report(all_metrics)
    report_path = os.path.join(MPNN_DIR, "comparison_report.md")
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"Report saved to: {report_path}")

    # Save all metrics as JSON
    metrics_path = os.path.join(MPNN_OUTPUT_DIR, "all_comparison_metrics.json")
    with open(metrics_path, 'w') as f:
        json.dump(all_metrics, f, indent=2, default=str)
    print(f"Metrics saved to: {metrics_path}")

    if failed_samples:
        print(f"\nFailed samples ({len(failed_samples)}):")
        for name, err in failed_samples:
            print(f"  {name}: {err}")

    print(f"\nDone! Processed {len(all_metrics)} samples successfully, {len(failed_samples)} failed.")
    print(f"Results saved to: {MPNN_OUTPUT_DIR}")
    print(f"Report saved to: {report_path}")


if __name__ == "__main__":
    main()
