# Example page
WaSH allows for arbitrary orderings of kernel functions in order to simulate a large variety of SPH scenarios. 

This page describes how one can use WaSH to create a basic water simulation.

# Initialising Common Parameters
The first step is to set some parameters for your simulation. 

```cpp
int main(int argc, char** argv) {
    wash::set_precision("double");
    wash::set_influence_radius(smoothingRadius);
    wash::set_max_iterations(1000);
}
```

# Outputs
The simulation results must be written *somewhere*. WaSH provides API calls for specifying what file to output to:
```cpp
    if (argc > 1) {
        // argv[1] = simulation name
        wash::set_simulation_name(argv[1]);
        if (argc > 2) {
            // argv[2] = output file name
            wash::set_output_file_name(argv[2]);
        } else {
            wash::set_output_file_name("ca");
        }
    } else {
        wash::set_simulation_name("serial_test");
    }
```

# Forces
A 'Force' in this case is not necessarily a force, but could be a particle's *property* for example temperature, that some other simulations may want to use for their calculations.
```cpp
    wash::add_force("nearDensity", 1);

    wash::add_force("predictedPosition", 2);
    wash::add_force("pressure", 2);
```

# Kernels
The main 
