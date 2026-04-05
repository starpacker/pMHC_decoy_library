import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    # Load data
    data_path = "data/GILGFVFTL_summary/Decoy_B/descriptor_solo_ranks.json"
    with open(data_path, 'r') as f:
        data = json.load(f)
    
    df = pd.DataFrame(data)
    
    # Extract rank columns
    rank_cols = [col for col in df.columns if col.startswith('rank_')]
    
    # Rename columns for better readability
    rename_map = {
        'rank_atchley_cosine': 'Atchley Cosine',
        'rank_sim_A': 'sim_A (MHC aligned)',
        'rank_sim_B': 'sim_B (Peptide aligned)',
        'rank_bf_corr': 'Boltz-2 Corr',
        'rank_rmsd_geo': 'RMSD Geo',
        'rank_esp_sim': 'ESP Sim',
        'rank_pesto_sim': 'PeSTo Sim'
    }
    
    df_ranks = df[rank_cols].rename(columns=rename_map)
    
    # Calculate Spearman correlation
    corr = df_ranks.corr(method='spearman')
    
    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    
    # 1. Heatmap of Spearman correlations
    sns.heatmap(corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1, center=0, 
                square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=axes[0])
    axes[0].set_title('Spearman Rank Correlation between Descriptors', fontsize=14)
    
    # 2. Parallel coordinates plot for top 10 by RMSD Geo
    top_10_rmsd = df.nsmallest(10, 'rank_rmsd_geo')
    top_10_ranks = top_10_rmsd[['sequence'] + rank_cols].rename(columns=rename_map)
    
    pd.plotting.parallel_coordinates(top_10_ranks, 'sequence', colormap='tab10', ax=axes[1], alpha=0.7, linewidth=2)
    axes[1].set_title('Rank Variations for Top 10 Candidates (by RMSD Geo)', fontsize=14)
    axes[1].set_ylabel('Rank (1 is best)')
    axes[1].invert_yaxis() # Invert y-axis so rank 1 is at the top
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.tight_layout()
    
    # Save to decoy_b folder
    output_path = "decoy_b/solo_ranking_visualization.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_path}")

if __name__ == "__main__":
    main()
