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
When the simulation starts, it will use whatever is defined in the `init` function. The contents of the `init` function may look like this if your goal is to spawn uniformly distributed particles:
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
This is the last part of our Main function. It describes the high-level behaviour of the simulation, and all that's left is the low-level specification of what should happen in each kernel.

Notice that the kernels are to be specified for *one* particle. WaSH will take care of the looping and parallelisation.

### VelocityUpdate Implementation
This kernel function helps us predict where the particle will be in the next timestep.
```cpp
void VelocityUpdate(wash::Particle& particle) {
    particle.set_vel(particle.get_vel() + ExternelForces(particle.get_pos(), particle.get_vel()) * deltaTime);
    // std::cout << "Particle velocity: " << particle.get_vel() << std::endl;

    const double predictionFactor = 1 / 120.0;
    // set predicted pos to the real position + some timestep of current vel
    particle.set_pos(particle.get_force_vector("position") + particle.get_vel() * predictionFactor);
    // std::cout << "Particle pred pos: " << particle.get_pos() << std::endl;
}
```

### CalculateDensity Implementation
Calculating density of particles is common across most, if not all, SPH simulations. Here's how we'll define it for our example:

```cpp
void CalculateDensity(wash::Particle& particle, const std::vector<wash::Particle>& neighbours) {
    // std::cout << "Running Custom Density Func" << std::endl;
    double density = 1.0;
    double nearDensity = 1.0;

    for (auto& neighbour : neighbours) {
        auto offset = neighbour.get_pos() - particle.get_pos();
        double dst = offset.magnitude();

        density += DensityKernel(dst, smoothingRadius);
        nearDensity += NearDensityKernel(dst, smoothingRadius);
    }

    particle.set_density(density);
    particle.set_force_scalar("nearDensity", nearDensity);
}
```

### force_kernel Implementation
Here, forces are applied to particles to simulate incompressible fluid with some viscosity. In this case, we want something resembling water.
```cpp
void force_kernel(wash::Particle& particle, const std::vector<wash::Particle>& neighbours) {
    CalculatePressureForce(particle, neighbours);
    CalculateViscosity(particle, neighbours);
}
```

```cpp
void CalculatePressureForce(wash::Particle& particle, const std::vector<wash::Particle>& neighbours) {
    double density = particle.get_density();
    double nearDensity = particle.get_force_scalar("nearDensity");
    double pressure = PressureFromDensity(density);
    double nearPressure = NearPressureFromDensity(nearDensity);
    wash::Vec2D pressureForce = wash::Vec2D({0.0, 0.0});

    wash::Vec2D pos = particle.get_pos();

    for (auto& neighbour : neighbours) {
        wash::Vec2D neighbourPos = neighbour.get_pos();
        wash::Vec2D offsetToNeighbour = neighbourPos - pos;
        double dst = offsetToNeighbour.magnitude();
        
        wash::Vec2D dirToNeighbour = dst > 0.0 ? offsetToNeighbour / dst : wash::Vec2D({0.0, 1.0});

        double neighbourDensity = neighbour.get_density();
        double neighbourNearDensity = neighbour.get_force_scalar("nearDensity");
        double neighbourPressure = PressureFromDensity(neighbourDensity);
        double neighbourNearPressure = NearPressureFromDensity(neighbourNearDensity);

        double sharedPressure = (pressure + neighbourPressure) * 0.5;
        double sharedNearPressure = (nearPressure + neighbourNearPressure) * 0.5;

        pressureForce += dirToNeighbour * DensityDerivative(dst, smoothingRadius) * sharedPressure / neighbourDensity;
        // std::cout << "w density p " << pressureForce << std::endl;
        pressureForce +=
            dirToNeighbour * NearDensityDerivative(dst, smoothingRadius) * sharedNearPressure / neighbourNearDensity;
        // std::cout << "w near density p " << pressureForce << std::endl;
    }

    wash::Vec2D acceleration = pressureForce / density;
    // std::cout << "PRESSURE FORCE p" << pressureForce << std::endl;

    particle.set_force_vector("pressure", pressureForce / density);
    particle.set_vel(particle.get_vel() + acceleration * deltaTime);
}
```

```cpp
void CalculateViscosity(wash::Particle& particle, const std::vector<wash::Particle>& neighbours) {
    wash::Vec2D pos = particle.get_pos();

    wash::Vec2D viscosityForce = wash::Vec2D { 0.0, 0.0 };
    wash::Vec2D velocity = particle.get_vel();

    for (auto& neighbour : neighbours) {
        wash::Vec2D neighbourPos = neighbour.get_pos();
        wash::Vec2D offsetToNeighbour = neighbourPos - pos;
        double dst = offsetToNeighbour.magnitude();
        
        wash::Vec2D neighbourVelocity = neighbour.get_vel();
        viscosityForce += (neighbourVelocity - velocity) * ViscosityKernel(dst, smoothingRadius);
    }

    particle.set_force_vector("viscosity", viscosityForce * viscosityStrength);
    particle.set_vel(particle.get_vel() + viscosityForce * viscosityStrength * deltaTime);
}
```


### UpdatePositions Implementation
This one's pretty self-explanatory. Using basic SUVAT to update the particles' positions.
```cpp
void UpdatePositions(wash::Particle& particle) {
    // particle.set_pos(particle.get_pos() + particle.get_vel() * deltaTime);
    particle.set_force_vector("position", particle.get_force_vector("position") + particle.get_vel() * deltaTime);
}
```

### HandleCollisions Implementation
The HandleCollisions kernel ensures particles do not exit the bounds of the simulation, and instead bounce back against a 'wall' (with a little damping)
```cpp
void HandleCollisions(wash::Particle& particle) {
    wash::Vec2D pos = particle.get_force_vector("position");
    wash::Vec2D vel = particle.get_vel();

    const wash::Vec2D halfSize = boundsSize * 0.5;
    wash::Vec2D edgeDst = halfSize - pos.abs();

    if (*(edgeDst[0]) <= 0) {
        *(pos[0]) = halfSize.at(0) * wash::sgn(pos.at(0));
        *(vel[0]) *= -1 * collisionDamping;
    }

    if (*(edgeDst[1]) <= 0) {
        *(pos[1]) = halfSize.at(1) * wash::sgn(pos.at(1));
        *(vel[1]) *= -1 * collisionDamping;
    }
    // do any obstacle collision here

    particle.set_force_vector("position", pos);
    particle.set_vel(vel);
}
```