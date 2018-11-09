1. run `mcnp6 i=infinite_medium_mcnp.i o=infinite_medium_mcnp.o tasks 8`
2. run `export DATABASE_PATH=path-to-db`
3. run `./infinite_medium.py --db_path=$DATABASE_PATH --sim_name="infinite_medium" --num_particles=1e8 --threads=8`
4. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --mcnp_file=infinite_medium_mcnp.o --mcnp_file_start=5282 --mcnp_file_end=6280 --current`
[H Infinite Medium](h_infinite_medium_current.png "H Infinite Medium")