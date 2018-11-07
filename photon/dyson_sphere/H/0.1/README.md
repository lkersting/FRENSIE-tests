1. run `mcnp6 i=dyson_sphere_mcnp.i o=dyson_sphere_mcnp.o tasks 8`
2. run `export DATABASE_PATH=path-to-db`
3. run `./dyson_sphere.py --db_path=$DATABASE_PATH --sim_name="dyson_sphere" --num_particles=1e9 --thread=8`
4. run `./dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --mcnp_file=dyson_sphere_mcnp.o --mcnp_file_start=4246 --mcnp_file_end=5244 --flux --emin=0.0995 --emax=0.1`