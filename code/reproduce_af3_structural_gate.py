#!/usr/bin/env python3
"""Reproduce the AF3 structural-interactome gate summary used in the manuscript."""
from pathlib import Path
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', default='data/af3', help='Directory containing AF3 TSV tables')
    args = parser.parse_args()
    d = Path(args.data_dir)
    edges = pd.read_csv(d / 'structural_edges_w2_calibrated_v2.tsv', sep='	')
    null = pd.read_csv(d / 'route_breaking_null_summary_w2_calibrated_v2_high.tsv', sep='	')
    w2c = pd.read_csv(d / 'w2c_seed_stability_by_pair.tsv', sep='	')
    high = edges[edges['af3_confidence_gate'].eq('high')].copy()
    same_route = high['route_a'].eq(high['route_b'])
    decoy_high = int(((edges['pair_class'] == 'unrelated_decoy') & (edges['af3_confidence_gate'] == 'high')).sum())
    obs_same_edges = int(same_route.sum())
    obs_edges = int(high.shape[0])
    obs_fraction = obs_same_edges / obs_edges
    null_dict = dict(zip(null['metric'], null['value']))
    q975 = float(null_dict['null_q975'])
    empirical_p_key = 'empirical_p_ge_observed'
    if empirical_p_key not in null_dict:
        empirical_p_key = 'empirical_p_right_tail'
    empirical_p = float(null_dict[empirical_p_key])
    stable_high = w2c[w2c['seed_stability_gate'].eq('stable_high')]
    bridge_stable = stable_high[stable_high['pair_class'].eq('reciprocal_cross_route')]
    decoy_stable = stable_high[stable_high['pair_class'].eq('unrelated_decoy')]
    print(f'W2 pair predictions: {edges.shape[0]}')
    print(f'W2 high-confidence edges: {obs_edges}')
    print(f'W2 same-route high-confidence edges: {obs_same_edges}')
    print(f'W2 same-route fraction: {obs_fraction:.6f}')
    print(f'W2 unrelated-decoy high-pass count: {decoy_high}')
    print(f'W2 route-label null q97.5: {q975:.6f}')
    print(f'W2 empirical right-tail p: {empirical_p:.6f}')
    print(f'W2C stable-high Fe-S bridge candidates: {bridge_stable.shape[0]}')
    print(f'W2C stable-high unrelated decoys: {decoy_stable.shape[0]}')
    assert edges.shape[0] == 654
    assert obs_edges == 122
    assert obs_same_edges == 119
    assert decoy_high == 0
    assert round(obs_fraction, 3) == 0.975
    assert round(q975, 3) == 0.631
    assert empirical_p <= 0.0001
    assert bridge_stable.shape[0] == 2
    assert decoy_stable.shape[0] == 0


if __name__ == '__main__':
    main()
