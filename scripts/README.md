# README


To produce many shuffles

```bash
cd src-tauri
mkdir -p shuffles
for seed in $(seq 1 100); do
  ../src-tauri/target/release/phenoblend-batch --phenopackets ../../../data/hpo/mgd --hpo ../../../data/hpo/hp.json --hpoa ../../../data/hpo/phenotype.hpoa --out shuffles/shuffle_$seed.tsv --n-sim 1000 --shuffle-seed $seed
done
```
