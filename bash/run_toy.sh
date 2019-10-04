#!/bin/bash
python3.7 -m PCA.PCA_runner -d 10 -p ~/federated_DP_PCA/results/simulation_toy/ -f ~/federated_DP_PCA/data/toydata/toy_data.tsv -k 5 -s -r 20 -c -i 0 -v -u -t -A -B -C -e 0.01,0.1,0.5,1 -g 0.01,0.1,1


#python3.7 -m PCA.PCA_runner -d 10 -p ~/federated_DP_PCA/results/simulation_breast/suball/ -f ~/federated_DP_PCA/data/tcga_breast/data/data_sub_2.tsv,~/federated_DP_PCA/data/tcga_breast/data/data_sub_3.tsv,~/federated_DP_PCA/data/tcga_breast/data/data_sub_4.tsv,~/federated_DP_PCA/data/tcga_breast/data/data_sub_1.tsv -k 2 -s -r 1 -c -i 0 -v -u -t -D -E



