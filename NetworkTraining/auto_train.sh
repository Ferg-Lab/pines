#!/bin/bash

start_val=1
end_val=8

for (( dim = $start_val; dim <=$end_val; dim++  ))
do
	sed -i "s/parser.add_argument(\"--latent_dim\", default=.*/parser.add_argument(\"--latent_dim\", default=${dim}, type=int)/g" mmd_vae_biased_3HL.py
	python3 earlyStop_train_vae.py
	bash rename.sh ${dim}d
done
