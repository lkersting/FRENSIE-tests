1. run `mcnp6 i=dyson_sphere_mcnp.i o=dyson_sphere_mcnp.o tasks 8`

2. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=4245 --mcnp_file_end=5244 --output_file="results/mcnp-mu=0.989.txt"`
3. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=364035 --mcnp_file_end=365034 --output_file="results/mcnp-mu=0.707.txt"`
4. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=512775 --mcnp_file_end=513774 --output_file="results/mcnp-mu=0.499.txt"`
5. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=659505 --mcnp_file_end=660504 --output_file="results/mcnp-mu=0.260.txt"`
6. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=808245 --mcnp_file_end=809244 --output_file="results/mcnp-mu=0.0.txt"`
7. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=956985 --mcnp_file_end=957984 --output_file="results/mcnp-mu=-0.260.txt"
8. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=1103715 --mcnp_file_end=1104714 --output_file="results/mcnp-mu=-0.499.txt"`
9. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=1252455 --mcnp_file_end=1253454 --output_file="results/mcnp-mu=-0.707.txt"`
10. run `./dyson_sphere_extract_mcnp_data.py --mcnp_file="dyson_sphere_mcnp.o" --mcnp_file_start=1612245 --mcnp_file_end=1613244 --output_file="results/mcnp-mu=-0.989.txt"`

11. enter dhp (decoupled half-profile) directory and follow README

12. run `/dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --output_file="results/dhp-mu=0.989.txt"`
13. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=2 --output_file="results/dhp-mu=0.707.txt"`
14. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=3 --output_file="results/dhp-mu=0.499.txt"`
15. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=4 --output_file="results/dhp-mu=0.260.txt"`
16. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=5 --output_file="results/dhp-mu=0.0.txt"`
17. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=6 --output_file="results/dhp-mu=-0.260.txt"`
18. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=7 --output_file="results/dhp-mu=-0.499.txt"`
19. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=8 --output_file="results/dhp-mu=-0.707.txt"`
20. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dhp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=9 --output_file="results/dhp-mu=-0.989.txt"`

21. enter dfp (decoupled full-profile) directory and follow README

22. run `/dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --output_file="results/dfp-mu=0.989.txt"`
23. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=2 --output_file="results/dfp-mu=0.707.txt"`
24. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=3 --output_file="results/dfp-mu=0.499.txt"`
25. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=4 --output_file="results/dfp-mu=0.260.txt"`
26. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=5 --output_file="results/dfp-mu=0.0.txt"`
27. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=6 --output_file="results/dfp-mu=-0.260.txt"`
28. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=7 --output_file="results/dfp-mu=-0.499.txt"`
29. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=8 --output_file="results/dfp-mu=-0.707.txt"`
30. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="dfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=9 --output_file="results/dfp-mu=-0.989.txt"`

31. enter cfp (coupled full-profile) directory and follow README

32. run `/dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=1 --output_file="results/cfp-mu=0.989.txt"`
33. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=2 --output_file="results/cfp-mu=0.707.txt"`
34. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=3 --output_file="results/cfp-mu=0.499.txt"`
35. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=4 --output_file="results/cfp-mu=0.260.txt"`
36. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=5 --output_file="results/cfp-mu=0.0.txt"`
37. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=6 --output_file="results/cfp-mu=-0.260.txt"`
38. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=7 --output_file="results/cfp-mu=-0.499.txt"`
39. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=8 --output_file="results/cfp-mu=-0.707.txt"`
40. run `./dyson_sphere_extract_frensie_data.py --rendezvous_file="cfp/dyson_sphere_rendezvous_10.xml" --estimator_id=1 --entity_id=9 --output_file="results/cfp-mu=-0.989.txt"`

