# StateSpaceSimulation
A header which includes objects used to integrate, sample and record state space systems.

### Example - simulate "Airy" function:
```cpp
// define Airy function
struct Airy {
    void operator()(const std::vector<double>& x, std::vector<double>& xd, const double t) {
        xd[0] = x[1];
        xd[1] = -t * x[0];
    }
};

// define simulation parameters
Models::Airy system;
std::vector<double> x{ { 1.0, 0.0 } };
double t{};
double dt{ 0.1 };
double t_end{ 10.0 };

// run simulation
Simulation::RK4 integrator;
Simulation::History recorder;
while (t < t_end) {
    recorder({ t, x[0], x[1] });
    integrator(system, x, t, dt);
}

// export to csv
recorder.csv("airy", { "t", "x0", "x1" }); // generate a file of comma separated values
```

drawing `airy.csv`:
![Image](https://github.com/user-attachments/assets/04a86f64-6d3f-4074-88ae-621c2a4e99ee)
