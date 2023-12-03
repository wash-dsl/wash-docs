# Example Simulation
WaSH allows for arbitrary orderings of kernel functions in order to simulate a large variety of SPH scenarios. 

This page describes how one can use WaSH to create a basic water simulation. A finished version of this example may be found in the WaSH source code as the `ca_fluid_sim` example.

# Constants
Some constants are defined before starting the simulation for easy parameter tweaking. Here's what constants should be defined at the top:
```cpp
constexpr wash::Vec2D spawnCentre { 3.35, 0.51 };
constexpr wash::Vec2D initialVelocity { 0.0, 0.0 };
constexpr wash::Vec2D spawnSize { 7.0, 7.0 };
constexpr wash::Vec2D boundsSize { 17.1, 9.3 };

constexpr double jitterStr = 0.025;
constexpr double numParticles = 4032;
constexpr double gravity = -12.0;

constexpr double deltaTime = TIME_DELTA(1, 3);
constexpr double collisionDamping = 0.95;
constexpr double smoothingRadius = 0.35;

constexpr double targetDensity = 55.0;
constexpr double pressureMultiplier = 500.0;
constexpr double nearPressureMultiplier = 18.0;
constexpr double viscosityStrength = 0.06;
```

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
int main(int argc, char** argv) {
    ...
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
}
```

# Forces
A 'Force' in this case is not necessarily a force, but could be a particle's *property* for example temperature, that some other simulations may want to use for their calculations.
```cpp
int main(int argc, char** argv) {
    ...
    wash::add_force("nearDensity", 1);

    wash::add_force("predictedPosition", 2);
    wash::add_force("pressure", 2);
}
```
These forces are now tracked for each particle in the simulation and may be used within kernel functions.

# Kernels
The kernels that describe the specific computations to take place in the simulation must be registered with WaSH.

## Initialisation Kernel
Your simulation must have a starting state for your particles. This is done like so:
```cpp
    wash::add_init_kernel(&init);
```
When the simulation starts, it will use whatever is defined in the `init` function. The contents of the `init` function may look like this:
```cpp
void init() {
    std::cout << "Calculated Time Step: " << deltaTime << std::endl;

    SpawnParticles(spawnSize, numParticles);
}
```

```cpp
void SpawnParticles(const wash::Vec2D spawnSizeVec, const size_t particleCount) {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re(42);

    double s_x = spawnSizeVec.at(0);
    double s_y = spawnSizeVec.at(1);
    
    int numX = (int)std::ceil( std::sqrt(
        s_x / s_y * particleCount + (s_x - s_y) * (s_x - s_y) / (4 * s_y * s_y)
    ) - (s_x - s_y) / (2 * s_y));

    int numY = (int)std::ceil( (double)particleCount / (double)numX );
    int i = 0;

    for (int y = 0; y < numY; y++) {
        for (int x = 0; x < numX; x++) {
            if (i >= particleCount) break;

            double tx = numX <= 1 ? 0.5 : x / (numX - 1.0);
            double ty = numY <= 1 ? 0.5 : y / (numY - 1.0);

            double angle = unif(re) * PI * 2.0;
            wash::Vec2D dir = wash::Vec2D({ std::cos(angle), std::sin(angle) });
            wash::Vec2D jitter = dir * jitterStr * (unif(re) - 0.5);
            wash::Vec2D pos = wash::Vec2D({ (tx - 0.5) * s_x, (ty - 0.5) * s_y }) + jitter + spawnCentre;

            // wash::Particle newp = wash::Particle();
            // newp.set_force_vector("position", newp.get_pos());
            // newp.set_vel(initialVelocity);
            // // VelocityUpdate(newp); // call here as the first initial call before density kernel
            // wash::add_par(newp);
            auto& p = wash::create_particle(0.0, 1.0, smoothingRadius, pos, initialVelocity);
            p.set_force_vector("position", pos);

            if (i < 5) {
                std::cout << "Particle " << i << " position " << pos << std::endl;
            }
            
            i++;
        }
    }
}
```

## Force and Update Kernels
These kernels define how particles change. The registration order of these kernels is order-sensitive, meaning they are run (at each iteration) in the order specified. Here's how we will order our force kernels:
```cpp
int main(int argc, char** argv) {
    ...
    wash::add_update_kernel(&VelocityUpdate);
    wash::add_force_kernel(&CalculateDensity);
    wash::add_force_kernel(&force_kernel);
  
    wash::add_update_kernel(&UpdatePositions);
    wash::add_update_kernel(&HandleCollisions);

    wash::start();
}
```
This is the last part of our Main function, and we must write what computations we want our kernels to perform.


