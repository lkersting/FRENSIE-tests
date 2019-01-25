1. run `export DATABASE_PATH=path-to-db`
2. run `./infinite_medium.py --db_path=$DATABASE_PATH --sim_name="infinite_medium" --num_particles=1e8 --threads=8`
3. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=5339 --mcnp_file_end=6337 --flux --forward`
4. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=3 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=6343 --mcnp_file_end=7341 --flux --forward`
5. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=6 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=7347 --mcnp_file_end=8345 --flux --forward`
[H Infinite Medium](h_infinite_medium_current.png "H Infinite Medium")