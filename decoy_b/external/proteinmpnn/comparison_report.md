# Tfold vs ProteinMPNN: CDR3 Sequence Design Comparison Report

**Generated:** 2026-02-25 14:01:40

**Number of samples:** 54

**ProteinMPNN model:** v_48_020 (vanilla)

**Temperature:** 0.1 (for sampling)

**Seed:** 42


## Summary Statistics

| Metric | Mean | Std | Min | Max |
|--------|------|-----|-----|-----|
| Agreement Rate | 0.2306 | 0.0908 | 0.0625 | 0.4286 |
| KL Divergence (Tfold||MPNN) | 1.4980 | 0.1974 | 1.1145 | 2.0907 |
| Cosine Similarity | 0.4712 | 0.0670 | 0.2924 | 0.6085 |
| Tfold Mean Entropy | 1.2930 | 0.1108 | 1.0613 | 1.5072 |
| MPNN Mean Entropy | 2.2025 | 0.1317 | 1.8538 | 2.3824 |
| Tfold Mean Max Prob | 0.5832 | 0.0383 | 0.5131 | 0.6540 |
| MPNN Mean Max Prob | 0.3548 | 0.0465 | 0.2749 | 0.4698 |


## Interpretation

- **Agreement Rate**: Fraction of CDR3 positions where Tfold and ProteinMPNN predict the same amino acid.
- **KL Divergence**: Measures how different the probability distributions are (lower = more similar).
- **Cosine Similarity**: Measures similarity of probability vectors (higher = more similar).
- **Entropy**: Measures uncertainty/diversity of predictions (higher = more uncertain).
- **Max Prob**: Confidence of the top prediction (higher = more confident).


## Per-Sample Results

| Sample | Seq Len | #Masked | Agreement | KL Div | Cos Sim | Tfold AAs | MPNN AAs |
|--------|---------|---------|-----------|--------|---------|-----------|----------|
| 6zkw_E_D_A_B_C_sample_1 | 227 | 25 | 0.160 | 1.408 | 0.480 | ASTSAFPATYTTYTAVVRFATANGT | AAYDNGAGAKATTVQVYDNANNKTR |
| 6zkw_E_D_A_B_C_sample_2 | 227 | 25 | 0.160 | 1.618 | 0.477 | ASSSAFQTTYYYGYAVTSGGTNKIT | LSLNAAGGGSGTPVQVNDNALLKPR |
| 6zkw_E_D_A_B_C_sample_3 | 227 | 25 | 0.320 | 1.366 | 0.508 | ASSPVFAGAYNQYTAVSSDAAAKYT | QASDAGAGSSSSTVAVSDAALRKTV |
| 7dzn_D_E_A_B_C_sample_1 | 224 | 23 | 0.130 | 2.091 | 0.292 | ASQGYAYKGYIAFGGGDAAANTF | ASSDLLNKPVKVSDLLLGGKKPK |
| 7dzn_D_E_A_B_C_sample_2 | 224 | 23 | 0.130 | 1.574 | 0.366 | ASGHRAKKAHIALHTAAGYNNGQ | SSGGGGSQPVKVGNGQLGGKLPS |
| 7dzn_D_E_A_B_C_sample_3 | 224 | 23 | 0.130 | 1.721 | 0.373 | ASGAAAGKYYIAFDFGAAGNANY | ASYDLGGQPQKVYGLLLGSSSPR |
| 7n2o_F_D_A_B_C_sample_1 | 228 | 28 | 0.179 | 1.315 | 0.535 | ASSLRGGGNNNQEAARRGGQGGGGNNYT | ASPSDAGGAAATTQVPALRGAAMGEPTV |
| 7n2o_F_D_A_B_C_sample_2 | 228 | 28 | 0.250 | 1.465 | 0.495 | ASSGYGGNGSKGYAVRQSGGGGGNDFYI | ASPSDAGGGATTTQVPYSWGSAMGSQTV |
| 7n2o_F_D_A_B_C_sample_3 | 228 | 28 | 0.321 | 1.235 | 0.548 | ASSSFGGGATKYYAVRQNGGGPGFANEI | ASPSDGGGGSTPTQVPYTLGGAMGDQPV |
| 7n2q_G_E_H_I_J_sample_1 | 221 | 22 | 0.273 | 1.336 | 0.511 | ASSAGGGNFAKQTAVRNGDNFS | ASPSAGGAGATPVQVPDALQPV |
| 7n2q_G_E_H_I_J_sample_2 | 221 | 22 | 0.227 | 1.311 | 0.548 | ASSRTGGFTANQTAVRSGNFVR | ASPSAGGSGSTTEVVPDAGQTV |
| 7n2q_G_E_H_I_J_sample_3 | 221 | 22 | 0.364 | 1.409 | 0.522 | ASSRGGGGFAKGTAVGANNFYF | ASESGGGGGASPTVAESGGLPV |
| 7n2r_F_D_A_B_C_sample_1 | 222 | 22 | 0.318 | 1.200 | 0.579 | ASSTFAGTNGNQTAVSTRGNVI | ASPSAGGAAGTPTVAPTGGTPV |
| 7n2r_F_D_A_B_C_sample_2 | 222 | 22 | 0.318 | 1.242 | 0.554 | ASSAGGGSSATYYAVRPQNNLI | ASPSAGGASSALTAAPVGGQLV |
| 7n2r_F_D_A_B_C_sample_3 | 222 | 22 | 0.318 | 1.239 | 0.529 | ASSPGGGSNAKTTASRFGNNFY | ASPPAGGAGSAPTSAPTGGLPV |
| 7n2s_F_D_A_B_C_sample_1 | 227 | 27 | 0.222 | 1.317 | 0.530 | ASSTQGTGYETLTAVDAGGGGGGGQLT | ASPSDLGGGSLLELVPLSLGAALEPPQ |
| 7n2s_F_D_A_B_C_sample_2 | 227 | 27 | 0.259 | 1.349 | 0.522 | ASSTATGGANTQEAVRAAGSFGNNFFT | ASPSDGGGAATLTLVPTPLGAAGGPPQ |
| 7n2s_F_D_A_B_C_sample_3 | 227 | 27 | 0.333 | 1.397 | 0.516 | ASSTGQGSASTQTAVRRAGGGGNNNLI | ASPSDLGGAATLEAVPQLWGAAAGPLQ |
| 7n5c_E_D_A_B_C_sample_1 | 220 | 14 | 0.143 | 1.811 | 0.378 | ASGLYDAQFGYDYY | AAGCCSGPVCGLSS |
| 7n5c_E_D_A_B_C_sample_2 | 220 | 14 | 0.143 | 1.778 | 0.384 | ASGFSNAQYNYNYT | AAGDGSGPVCGLSS |
| 7n5c_E_D_A_B_C_sample_3 | 220 | 14 | 0.071 | 1.813 | 0.365 | ASGVYYAQFNNDTY | SAGCCSGPVSGPSS |
| 7na5_E_D_A_B_C_sample_1 | 223 | 22 | 0.136 | 1.510 | 0.432 | AADRNFGSQYVQTAVGSLAKFI | VAYPDGGGGKGPVAAYDDGRPE |
| 7na5_E_D_A_B_C_sample_2 | 223 | 22 | 0.136 | 1.590 | 0.423 | AWDGQALGPYKWFAAQTYAQLD | VAVPGGAGGKLPKAAVGVGPPQ |
| 7na5_E_D_A_B_C_sample_3 | 223 | 22 | 0.273 | 1.477 | 0.432 | AWSSNPGGGNGYFAVGFGAGLI | AAYPGGGGGKGIVAQYGVGPIV |
| 7nme_E_D_A_B_C_sample_1 | 221 | 20 | 0.200 | 1.489 | 0.466 | ASGASYTQFATTNYNNNKYR | ASGDLLGLQAAGSGTAGFPV |
| 7nme_E_D_A_B_C_sample_2 | 221 | 20 | 0.350 | 1.442 | 0.501 | ASGRHGKMFAATASNSSKYI | ASGDLGGPQAAGSSSAGGPV |
| 7nme_E_D_A_B_C_sample_3 | 221 | 20 | 0.350 | 1.370 | 0.496 | ASGRADKQFAADDFNTGGTI | ASGDLLGLQAAGSGVAGGPV |
| 7ow6_E_D_A_B_C_sample_1 | 232 | 28 | 0.214 | 1.537 | 0.424 | ASRVVGGAGGYTYTALSGVQGGGRNKLI | AAAPVVVAADEGIQAIAPVVVAGGLQIV |
| 7ow6_E_D_A_B_C_sample_2 | 232 | 28 | 0.214 | 1.551 | 0.420 | ASSVVNPGSGFSYFAFSAFLGGGGYTLI | AASGAAAAADEALQAISGPAAAGSLRLV |
| 7ow6_E_D_A_B_C_sample_3 | 232 | 28 | 0.143 | 1.788 | 0.337 | ASGVDRDGGAGTYHALSAFAQGAGNKFI | LAASAAAGAAEQLQVIAGAAAAADLLLV |
| 7pbe_E_D_A_B_C_sample_1 | 222 | 21 | 0.286 | 1.512 | 0.463 | ASGRDGGGKYFVGYYGGGNFI | ASESGSGGPPVAVPDGSGFPS |
| 7pbe_E_D_A_B_C_sample_2 | 222 | 21 | 0.286 | 1.545 | 0.471 | ASGALGGGAYFVEGAYDGGMT | ASSSGSGGGPVLVDSGDGFPE |
| 7pbe_E_D_A_B_C_sample_3 | 222 | 21 | 0.286 | 1.536 | 0.425 | ASGLEGGGKYFVGGAGDNNMI | ASEGGSGGPPVAVPDGDGPPE |
| 7phr_B_A_H_L_P_sample_1 | 221 | 20 | 0.250 | 1.676 | 0.506 | AGGLYGGSGQYATGGAGKLD | AAASGGGDAIVSAAGGGPIE |
| 7phr_B_A_H_L_P_sample_2 | 221 | 20 | 0.150 | 1.691 | 0.476 | AGGGYANGGQYATYDAAKLD | AAGSGGGDGPVLAGGGGGPE |
| 7phr_B_A_H_L_P_sample_3 | 221 | 20 | 0.200 | 1.598 | 0.478 | ASGGAYNGGQFATYSGGKFI | AAGSPAADAPVVAGAGGPPV |
| 7qpj_B_A_C_D_E_sample_1 | 221 | 20 | 0.400 | 1.708 | 0.438 | ASSLLGYYFAVGTGGGNNLT | AAGGSGGPQAVGNGGGDPPK |
| 7qpj_B_A_C_D_E_sample_2 | 221 | 20 | 0.300 | 1.770 | 0.455 | ASRFADGYYAVTVDGTFKLQ | ASGGGGGIEAVGNGGGDGIK |
| 7qpj_B_A_C_D_E_sample_3 | 221 | 20 | 0.400 | 1.505 | 0.499 | ASGDAGGGYAVGGDAWAGFQ | AAVGAGGPVAVVNGAGDGPK |
| 7r80_B_A_C_D_E_sample_1 | 226 | 24 | 0.167 | 1.485 | 0.474 | ASSPFGGGFYYGTARVYFGAAKYT | AAWSGSSGSTPITAVPSSSGLGIR |
| 7r80_B_A_C_D_E_sample_2 | 226 | 24 | 0.208 | 1.602 | 0.402 | ASSYLLLPGTNATAVGAFYTAQFT | AAWDGAAGSTAVTAVWSGAGLGVR |
| 7r80_B_A_C_D_E_sample_3 | 226 | 24 | 0.333 | 1.479 | 0.517 | ASSYGAGGGTTYTAVRAYGSAGLT | AAWGAAGGALGITAVWGSAGGGIS |
| 7rk7_E_D_A_B_C_sample_1 | 232 | 32 | 0.094 | 1.565 | 0.432 | AIGSRAGIGGITYVTIQHLAGGSRFGANYNGI | VSVGNNNNNGGTSLGQIRVVVGNNAATTTTIT |
| 7rk7_E_D_A_B_C_sample_2 | 232 | 32 | 0.156 | 1.463 | 0.435 | AIGSRHSGIGDTEDGTQFGAGGRDAGAGSGFI | VSLGGGGGGGAKGLGEPRVVLGGGTGGSGLPH |
| 7rk7_E_D_A_B_C_sample_3 | 232 | 32 | 0.062 | 1.638 | 0.432 | AISERAGYSGGETVNSQTLLPGDDDLAFAKQI | VSLGGGGGGSAPPFGEPRAVLGGGAATGGLPV |
| 7t2b_O_N_K_L_M_sample_1 | 224 | 23 | 0.261 | 1.518 | 0.470 | ASPGDGEFYHATGAGAGGNYGLI | AAGGGGGRLVLAPSGSSGGLPLS |
| 7t2b_O_N_K_L_M_sample_2 | 224 | 23 | 0.087 | 1.696 | 0.392 | ASGRGGEAGEATRAAAAGSNGLN | AAPGSSGRPVKAPTGSSGGLLPS |
| 7t2b_O_N_K_L_M_sample_3 | 224 | 23 | 0.130 | 1.756 | 0.385 | ASGQAAREGYATGAAGTYTFKLT | KAGGASGREVKAPTAPAGGSLES |
| 7z50_E_H_A_B_T_sample_1 | 223 | 21 | 0.190 | 1.114 | 0.609 | ASSAAGGGTYYAASAFSVSII | GSASLLGLALEAAALAGARLH |
| 7z50_E_H_A_B_T_sample_2 | 223 | 21 | 0.190 | 1.313 | 0.543 | ASGGADSNTWFAASAAWGALI | ASPLGELLLPEAAPLGGLKPV |
| 7z50_E_H_A_B_T_sample_3 | 223 | 21 | 0.190 | 1.133 | 0.550 | ASGASHRNTYTAATAAYASFT | ASPSGGSLGPEAAPLGGGKPK |
| 8cx4_F_D_A_B_C_sample_1 | 230 | 28 | 0.286 | 1.257 | 0.569 | AGSTIGGTNNKLYAVRLSGGSGGANFLT | ASPGDGGAGGQIVAVPGLLGAGMAEPIQ |
| 8cx4_F_D_A_B_C_sample_2 | 230 | 28 | 0.321 | 1.324 | 0.523 | ASSSTGGGANPQTAVRIQGGGGGANAFT | ASPGGGGGESGIVAVPGLGGAAMGEPIQ |
| 8cx4_F_D_A_B_C_sample_3 | 230 | 28 | 0.429 | 1.255 | 0.555 | ASSAAGGGQNTYTAVRQYGTGGGGFGLT | ASPSDGGGGSTLSAVPTLGGAGMGEPLQ |


## Detailed Per-Position Comparison (First 5 Samples)


### 6zkw_E_D_A_B_C_sample_1

- Tfold CDR3: `ASTSAFPATYTTYTAVVRFATANGT`
- MPNN  CDR3: `AAYDNGAGAKATTVQVYDNANNKTR`
- Agreement: 0.160


**Chain B:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 91 | A | R | ✗ | 3.206 | 0.113 | 1.000 | 0.217 |
| 92 | S | A | ✗ | 1.285 | 0.404 | 0.911 | 0.501 |
| 93 | T | Y | ✗ | 2.014 | 0.190 | 0.527 | 0.224 |
| 94 | S | P | ✗ | 0.700 | 0.637 | 0.265 | 0.099 |
| 95 | A | R | ✗ | 1.713 | 0.385 | 0.829 | 0.095 |
| 96 | F | A | ✗ | 1.112 | 0.458 | 0.305 | 0.410 |
| 97 | P | A | ✗ | 1.833 | 0.308 | 0.457 | 0.499 |
| 98 | A | G | ✗ | 0.971 | 0.808 | 0.381 | 0.734 |
| 99 | T | S | ✗ | 0.641 | 0.648 | 0.286 | 0.169 |
| 100 | Y | T | ✗ | 1.059 | 0.574 | 0.258 | 0.116 |
| 101 | T | A | ✗ | 1.316 | 0.364 | 0.404 | 0.147 |
| 102 | T | T | ✓ | 0.780 | 0.683 | 0.409 | 0.155 |
| 103 | Y | T | ✗ | 2.149 | 0.255 | 0.483 | 0.177 |
| 104 | T | Q | ✗ | 1.529 | 0.427 | 0.535 | 0.132 |

**Chain A:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 90 | A | Q | ✗ | 1.734 | 0.510 | 0.999 | 0.189 |
| 91 | V | V | ✓ | 0.639 | 0.936 | 0.972 | 0.480 |
| 92 | V | T | ✗ | 0.907 | 0.746 | 0.268 | 0.439 |
| 93 | R | S | ✗ | 1.945 | 0.263 | 0.319 | 0.374 |
| 94 | F | A | ✗ | 1.549 | 0.352 | 0.309 | 0.181 |
| 95 | A | E | ✗ | 1.575 | 0.322 | 0.405 | 0.244 |
| 96 | T | D | ✗ | 1.263 | 0.360 | 0.460 | 0.164 |
| 97 | A | R | ✗ | 1.343 | 0.405 | 0.654 | 0.129 |
| 98 | N | K | ✗ | 0.983 | 0.870 | 0.418 | 0.600 |
| 99 | G | Y | ✗ | 0.845 | 0.716 | 0.427 | 0.200 |
| 100 | T | R | ✗ | 2.112 | 0.267 | 0.903 | 0.198 |

### 6zkw_E_D_A_B_C_sample_2

- Tfold CDR3: `ASSSAFQTTYYYGYAVTSGGTNKIT`
- MPNN  CDR3: `LSLNAAGGGSGTPVQVNDNALLKPR`
- Agreement: 0.160


**Chain B:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 91 | A | Y | ✗ | 4.808 | 0.019 | 1.000 | 0.348 |
| 92 | S | A | ✗ | 1.286 | 0.428 | 0.911 | 0.443 |
| 93 | S | T | ✗ | 0.892 | 0.792 | 0.527 | 0.295 |
| 94 | S | P | ✗ | 0.528 | 0.776 | 0.265 | 0.145 |
| 95 | A | G | ✗ | 2.081 | 0.196 | 0.829 | 0.311 |
| 96 | F | G | ✗ | 1.466 | 0.718 | 0.305 | 0.830 |
| 97 | Q | G | ✗ | 1.780 | 0.367 | 0.457 | 0.663 |
| 98 | T | G | ✗ | 0.846 | 0.812 | 0.381 | 0.664 |
| 99 | T | S | ✗ | 0.718 | 0.614 | 0.286 | 0.156 |
| 100 | Y | S | ✗ | 1.287 | 0.404 | 0.258 | 0.138 |
| 101 | Y | G | ✗ | 1.488 | 0.306 | 0.404 | 0.162 |
| 102 | Y | T | ✗ | 0.770 | 0.632 | 0.409 | 0.123 |
| 103 | G | I | ✗ | 4.320 | 0.031 | 0.483 | 0.925 |
| 104 | Y | V | ✗ | 1.479 | 0.412 | 0.535 | 0.128 |

**Chain A:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 90 | A | A | ✓ | 0.833 | 0.833 | 0.999 | 0.432 |
| 91 | V | V | ✓ | 0.353 | 0.982 | 0.972 | 0.653 |
| 92 | T | S | ✗ | 1.579 | 0.537 | 0.268 | 0.803 |
| 93 | S | D | ✗ | 1.911 | 0.176 | 0.319 | 0.503 |
| 94 | G | S | ✗ | 1.397 | 0.369 | 0.309 | 0.171 |
| 95 | G | A | ✗ | 1.668 | 0.394 | 0.405 | 0.199 |
| 96 | T | L | ✗ | 1.481 | 0.208 | 0.460 | 0.318 |
| 97 | N | A | ✗ | 0.943 | 0.641 | 0.654 | 0.168 |
| 98 | K | K | ✓ | 0.670 | 0.842 | 0.418 | 0.297 |
| 99 | I | G | ✗ | 4.065 | 0.013 | 0.427 | 0.899 |
| 100 | T | R | ✗ | 1.794 | 0.420 | 0.903 | 0.126 |

### 6zkw_E_D_A_B_C_sample_3

- Tfold CDR3: `ASSPVFAGAYNQYTAVSSDAAAKYT`
- MPNN  CDR3: `QASDAGAGSSSSTVAVSDAALRKTV`
- Agreement: 0.320


**Chain B:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 91 | A | L | ✗ | 3.108 | 0.106 | 1.000 | 0.356 |
| 92 | S | A | ✗ | 1.045 | 0.569 | 0.911 | 0.388 |
| 93 | S | S | ✓ | 1.266 | 0.502 | 0.527 | 0.226 |
| 94 | P | P | ✓ | 0.584 | 0.773 | 0.265 | 0.188 |
| 95 | V | S | ✗ | 1.511 | 0.455 | 0.829 | 0.133 |
| 96 | F | A | ✗ | 2.827 | 0.099 | 0.305 | 0.907 |
| 97 | A | A | ✓ | 2.043 | 0.261 | 0.457 | 0.613 |
| 98 | G | G | ✓ | 0.721 | 0.834 | 0.381 | 0.596 |
| 99 | A | S | ✗ | 0.769 | 0.578 | 0.286 | 0.146 |
| 100 | Y | A | ✗ | 1.185 | 0.542 | 0.258 | 0.277 |
| 101 | N | G | ✗ | 1.275 | 0.380 | 0.404 | 0.118 |
| 102 | Q | T | ✗ | 0.805 | 0.645 | 0.409 | 0.139 |
| 103 | Y | Y | ✓ | 0.957 | 0.817 | 0.483 | 0.603 |
| 104 | T | Q | ✗ | 1.608 | 0.367 | 0.535 | 0.123 |

**Chain A:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 90 | A | L | ✗ | 2.299 | 0.300 | 0.999 | 0.197 |
| 91 | V | V | ✓ | 0.381 | 0.980 | 0.972 | 0.632 |
| 92 | S | S | ✓ | 0.625 | 0.742 | 0.268 | 0.303 |
| 93 | S | D | ✗ | 1.250 | 0.392 | 0.319 | 0.220 |
| 94 | D | V | ✗ | 1.260 | 0.436 | 0.309 | 0.121 |
| 95 | A | A | ✓ | 1.618 | 0.478 | 0.405 | 0.389 |
| 96 | A | L | ✗ | 1.609 | 0.218 | 0.460 | 0.270 |
| 97 | A | V | ✗ | 1.546 | 0.289 | 0.654 | 0.200 |
| 98 | K | K | ✓ | 0.950 | 0.784 | 0.418 | 0.321 |
| 99 | Y | Y | ✓ | 1.124 | 0.752 | 0.427 | 0.722 |
| 100 | T | V | ✗ | 1.791 | 0.409 | 0.903 | 0.157 |

### 7dzn_D_E_A_B_C_sample_1

- Tfold CDR3: `ASQGYAYKGYIAFGGGDAAANTF`
- MPNN  CDR3: `ASSDLLNKPVKVSDLLLGGKKPK`
- Agreement: 0.130


**Chain B:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 92 | A | S | ✗ | 1.734 | 0.486 | 0.999 | 0.260 |
| 93 | S | S | ✓ | 0.376 | 0.972 | 0.872 | 0.734 |
| 94 | Q | F | ✗ | 3.823 | 0.013 | 0.539 | 0.897 |
| 95 | G | G | ✓ | 1.351 | 0.398 | 0.282 | 0.253 |
| 96 | Y | G | ✗ | 2.942 | 0.150 | 0.541 | 0.904 |
| 97 | A | G | ✗ | 2.172 | 0.227 | 0.778 | 0.140 |
| 98 | Y | G | ✗ | 0.370 | 0.783 | 0.239 | 0.146 |
| 99 | K | P | ✗ | 1.778 | 0.277 | 0.782 | 0.217 |
| 100 | G | P | ✗ | 2.946 | 0.025 | 0.209 | 0.824 |
| 101 | Y | V | ✗ | 1.629 | 0.258 | 0.495 | 0.248 |

**Chain A:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 88 | I | V | ✗ | 2.561 | 0.144 | 0.932 | 0.362 |
| 89 | A | V | ✗ | 4.090 | 0.054 | 0.936 | 0.538 |
| 90 | F | Q | ✗ | 3.516 | 0.014 | 0.455 | 0.890 |
| 91 | G | G | ✓ | 2.093 | 0.244 | 0.533 | 0.768 |
| 92 | G | Y | ✗ | 3.032 | 0.035 | 0.293 | 0.882 |
| 93 | G | Y | ✗ | 2.397 | 0.090 | 0.552 | 0.746 |
| 94 | D | L | ✗ | 2.867 | 0.087 | 0.709 | 0.296 |
| 95 | A | G | ✗ | 0.736 | 0.519 | 0.632 | 0.515 |
| 96 | A | S | ✗ | 1.037 | 0.541 | 0.421 | 0.129 |
| 97 | A | K | ✗ | 1.753 | 0.276 | 0.600 | 0.144 |
| 98 | N | S | ✗ | 1.064 | 0.485 | 0.309 | 0.132 |
| 99 | T | G | ✗ | 1.505 | 0.474 | 0.273 | 0.660 |
| 100 | F | A | ✗ | 2.313 | 0.175 | 0.466 | 0.123 |

### 7dzn_D_E_A_B_C_sample_2

- Tfold CDR3: `ASGHRAKKAHIALHTAAGYNNGQ`
- MPNN  CDR3: `SSGGGGSQPVKVGNGQLGGKLPS`
- Agreement: 0.130


**Chain B:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 92 | A | A | ✓ | 0.702 | 0.850 | 0.999 | 0.494 |
| 93 | S | S | ✓ | 0.388 | 0.962 | 0.872 | 0.702 |
| 94 | G | L | ✗ | 1.849 | 0.152 | 0.539 | 0.436 |
| 95 | H | D | ✗ | 0.939 | 0.502 | 0.282 | 0.235 |
| 96 | R | T | ✗ | 0.844 | 0.452 | 0.541 | 0.298 |
| 97 | A | G | ✗ | 2.028 | 0.250 | 0.778 | 0.159 |
| 98 | K | E | ✗ | 0.471 | 0.722 | 0.239 | 0.133 |
| 99 | K | P | ✗ | 1.863 | 0.241 | 0.782 | 0.246 |
| 100 | A | P | ✗ | 1.411 | 0.245 | 0.209 | 0.273 |
| 101 | H | V | ✗ | 1.741 | 0.255 | 0.495 | 0.206 |

**Chain A:**

| Pos | Tfold AA | MPNN AA | Match | KL Div | Cos Sim | Tfold Conf | MPNN Conf |
|-----|----------|---------|-------|--------|---------|------------|-----------|
| 88 | I | K | ✗ | 3.358 | 0.071 | 0.932 | 0.283 |
| 89 | A | V | ✗ | 3.678 | 0.065 | 0.936 | 0.480 |
| 90 | L | G | ✗ | 1.997 | 0.155 | 0.455 | 0.436 |
| 91 | H | N | ✗ | 1.084 | 0.528 | 0.533 | 0.178 |
| 92 | T | R | ✗ | 2.154 | 0.095 | 0.293 | 0.657 |
| 93 | A | R | ✗ | 1.772 | 0.213 | 0.552 | 0.240 |
| 94 | A | L | ✗ | 2.146 | 0.203 | 0.709 | 0.183 |
| 95 | G | G | ✓ | 0.888 | 0.494 | 0.632 | 0.432 |
| 96 | Y | G | ✗ | 0.851 | 0.629 | 0.421 | 0.157 |
| 97 | N | K | ✗ | 1.647 | 0.298 | 0.600 | 0.126 |
| 98 | N | S | ✗ | 1.123 | 0.464 | 0.309 | 0.120 |
| 99 | G | A | ✗ | 1.309 | 0.365 | 0.273 | 0.156 |
| 100 | Q | S | ✗ | 1.959 | 0.216 | 0.466 | 0.122 |
