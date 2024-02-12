import argparse
import scanpy as sc
import numpy as np
import pandas as pd

import pyVIA.core as via
import matplotlib
matplotlib.use('Agg')
  
  
if __name__ == '__main__':
  p = argparse.ArgumentParser()
  p.add_argument('--data_file', type=str, default="BEYOND/data/BEYOND.DLPFC.h5ad")
  p.add_argument('--output_path', type=str, default="all/data/via.pickle")
  p.add_argument('--k', help='Number of neighbors for neighboring graph over which clustering is performed', type=int, required=True)
  p.add_argument('--big', type=float, nargs = "+", required=True)
  p.add_argument('--small', type=int, nargs = "+", required=True)
  p.add_argument('--root_clusters', type=str, nargs="+", required=True)
  p.add_argument('--exclude_clusters', type=str, nargs="*", required=False, default=None)
  args = p.parse_args()

  print('----------------------------------------------------------')
  print(f"Running VIA with the following arguments:\n\tData file: {args.data_file}\n\tOutput file: {args.output_path}\n\tk: {args.k}"\
  f"\n\tToo big fraction: {args.big}\n\tSmall cluster threshold: {args.small}\n\tRoot clusters: {args.root_clusters}\n\tExcluded clusters: {args.exclude_clusters}")
  print('----------------------------------------------------------\n')

  data = sc.read_h5ad(args.data_file)
  full_obs = pd.DataFrame(index=data.obs_names.values)
  
  if args.exclude_clusters is not None:
    data = data[~data.obs["clusters"].isin(args.exclude_clusters)]

  mass = np.array(data[data.obs["clusters"].isin(args.root_clusters)].X)
  dist = np.linalg.norm(mass - np.mean(mass, axis=0), axis=1)
  root = np.where(dist == dist.min())[0]

  v0 = via.VIA(data.X, true_label = data.obs["clusters"], root_user=root,
               knn=args.k, too_big_factor=args.big[0], small_pop=args.small[0],
               is_coarse=True, preserve_disconnected=True, n_iter_leiden=50, num_threads=1)
  v0.run_VIA()

  ts = via._get_loc_terminal_states(v0, data.X)
  v  = via.VIA(data.X, true_label = data.obs["clusters"], root_user=root,
               knn=args.k, too_big_factor=args.big[-1], small_pop=args.small[-1],
               num_mcmc_simulations = 2500,
               super_cluster_labels=v0.labels,
               super_terminal_cells=ts,
               super_node_degree_list=v0.node_degree_list,
               super_terminal_clusters=v0.terminal_clusters,
               full_neighbor_array=v0.full_neighbor_array,
               full_distance_array=v0.full_distance_array,
               ig_full_graph=v0.ig_full_graph,
               csr_array_locally_pruned=v0.csr_array_locally_pruned,
               is_coarse=False, preserve_disconnected=True, num_threads=1, n_iter_leiden=50)
  v.run_VIA()
  ts = via._get_loc_terminal_states(v, data.X)
  
  
  in_full = lambda df: full_obs.merge(pd.DataFrame(df, index=data.obs_names.values), how="left", left_index=True, right_index=True)

  import pickle
  with open(args.output_path, 'wb') as handle:
    pickle.dump({"idents": data.obs_names.values, 
                 "user.root":data.obs_names[root], 
                 "pseudotime": in_full(v.single_cell_pt_markov),
                 "trajectories.scaled": in_full(v.single_cell_bp),
                 "trajectories.normalized": in_full(v.single_cell_bp_rownormed),
                 "terminals": ts}, 
      handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved VIA trajectories output to {args.output_path}")
