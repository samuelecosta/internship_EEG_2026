Working setup:

```
main_dir/
├── ERCEM321_Atreyu/
│   ├── build/
│   │   ├── bin/Release/
│   │   └── lib/
│   └── Examples/
│       └── matlab_shared_lib_ThreeLayers/
│           ├── aux/
│           ├── mesh/
│           └── matlab_shared_lib_ThreeLayers.h
│
└── my_dir/
    ├── scripts/
    │   └── matlab/
    └── datasets/
        └── tvb_default/
            ├── h5_files/
            └── msh_files/
```

The TVB python code has to be used inside TVB jupyter notebook

Workflow:

Generate the surfaces with surface_utils.m -> Generate the matrices with forward_problem_const.m and three_layer_gen.m -> Compare the two with G_comparison_scalp.m -> See inverse algs performance with inverse_problem_validation.m -> Create a tvb simulation with surface_simulation.py and use the generated data to complete the validation pipeline in tvb_validation.m
