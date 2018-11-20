1. run `export DATABASE_PATH=path-to-db`
2. run `./dyson_sphere.py --db_path=$DATABASE_PATH --sim_name="dyson_sphere" --num_particles=1e9 --thread=8`

# Phi = 0.15
3. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=1770 --mcnp_file_end=2768 --emin=0.0995 --emax=0.1 --top_ymax=3e-8 --bottom_ymin=0.60 --bottom_ymax=1.40 --legend_xpos=1.01 --legend_ypos=1.05`

4. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=2 --entity_id=1 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=40050 --mcnp_file_end=50030 --emin=0.09955 --emax=0.1 --top_ymax=4e-8 --bottom_ymin=0.60 --bottom_ymax=1.40 --legend_xpos=1.01 --legend_ypos=1.05`

# Phi = Pi/12
5. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=2 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=4785 --mcnp_file_end=5783 --emin=0.0988 --emax=0.1 --top_ymax=2e-8 --bottom_ymin=0.9 --bottom_ymax=1.1 --legend_xpos=1.01 --legend_ypos=1.05`

6. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=2 --entity_id=2 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=50036 --mcnp_file_end=60016 --emin=0.0988 --emax=0.1 --top_ymax=2.5e-8 --bottom_ymin=0.6 --bottom_ymax=1.4 --legend_xpos=1.01 --legend_ypos=1.05`

# Phi = Pi/6
7. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=3 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=7800 --mcnp_file_end=8798 --emin=0.096 --emax=0.099 --top_ymax=1e-8 --bottom_ymin=0.6 --bottom_ymax=1.4 --legend_xpos=1.01 --legend_ypos=1.05`

8. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=2 --entity_id=3 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=60022 --mcnp_file_end=70002 --emin=0.096 --emax=0.099 --top_ymax=1e-8 --bottom_ymin=0.6 --bottom_ymax=1.4 --legend_xpos=1.01 --legend_ypos=1.05`

# Phi = Pi/4
9. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=4 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=10815 --mcnp_file_end=11813 --emin=0.093 --emax=0.096 --top_ymax=6e-9 --bottom_ymin=0.6 --bottom_ymax=1.4 --legend_xpos=1.01 --legend_ypos=1.05`

10. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=2 --entity_id=4 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=70008 --mcnp_file_end=79988 --emin=0.093 --emax=0.096 --top_ymax=6e-9 --bottom_ymin=0.6 --bottom_ymax=1.4 --legend_xpos=1.01 --legend_ypos=1.05`

5. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=3 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=512776 --mcnp_file_end=513774 --emin=0.089 --emax=0.093 --top_ymax=4e-9 --bottom_ymin=0.8 --bottom_ymax=1.2 --legend_xpos=1.00 --legend_ypos=1.00`
6. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=4 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=659506 --mcnp_file_end=660504 --emin=0.085 --emax=0.090 --top_ymax=2.5e-9 --bottom_ymin=0.8 --bottom_ymax=1.2 --legend_xpos=1.00 --legend_ypos=1.02`
7. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=5 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=808246 --mcnp_file_end=809244 --emin=0.081 --emax=0.087 --top_ymax=2e-9 --bottom_ymin=0.8 --bottom_ymax=1.2 --legend_xpos=1.00 --legend_ypos=1.02`
8. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=6 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=956986 --mcnp_file_end=957984 --emin=0.0785 --emax=0.082 --top_ymax=2e-9 --bottom_ymin=0.8 --bottom_ymax=1.2 --legend_xpos=1.00 --legend_ypos=1.02`
9. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=7 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=1103716 --mcnp_file_end=1104714 --emin=0.075 --emax=0.080 --top_ymax=2e-9 --bottom_ymin=0.8 --bottom_ymax=1.2 --legend_xpos=1.00 --legend_ypos=1.02`
10. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=8 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=1252456 --mcnp_file_end=1253454 --emin=0.072 --emax=0.078 --top_ymax=2.5e-9 --bottom_ymin=0.8 --bottom_ymax=1.2 --legend_xpos=1.00 --legend_ypos=1.02`
11. run `../dyson_sphere-plot.py --rendezvous_file="dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=9 --mcnp_file=../dyson_sphere_mcnp.o --mcnp_file_start=1612246 --mcnp_file_end=1613244 --emin=0.070 --emax=0.074 --top_ymax=3e-9 --bottom_ymin=0.8 --bottom_ymax=1.2 --legend_xpos=1.00 --legend_ypos=1.02`

# ../dyson_sphere_restart.py --db_path=$DATABASE_PATH --rendezvous_file="dyson_sphere_rendezvous_10.xml" --num_particles=1e9 --threads=8