1. run `export DATABASE_PATH=path-to-db`
2. run `./infinite_medium.py --db_path=$DATABASE_PATH --sim_name="infinite_medium" --num_particles=1e8 --threads=8`
3. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=5337 --mcnp_file_end=6335 --flux --forward`
4. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=3 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=6341 --mcnp_file_end=7339 --flux --forward`
5. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=6 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=7345 --mcnp_file_end=8343 --flux --forward`
5. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=9 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=8349 --mcnp_file_end=9347 --flux --forward`
6. run `./infinite_medium-plot.py --rendezvous_file="infinite_medium_rendezvous_10.xml" --estimator_id=1 --entity_id=12 --mcnp_file=../infinite_medium_mcnp.o --mcnp_file_start=9353 --mcnp_file_end=10351 --flux --forward`
[H Infinite Medium](h_infinite_medium_current.png "H Infinite Medium")