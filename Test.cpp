#include "SimulinkCPP.hpp"

/**
* models to be simulated
**/
namespace Models {
    /**
    * \brief Airy function (https://en.wikipedia.org/wiki/Airy_function)
    **/
    struct Airy {
        void operator()(const std::vector<double>& x, std::vector<double>& xd, const double t) {
            xd[0] = x[1];
            xd[1] = -t * x[0];
        }
    };

    /**
    * \brief Lorenz attractor (https://en.wikipedia.org/wiki/Lorenz_system)
    **/
    struct Lorenz {
        static constexpr double sigma{ 10.0 };
        static constexpr double R{ 28.0 };
        static constexpr double b{ 8.0 / 3.0 };

        void operator()(const std::vector<double>& x, std::vector<double>& xd, const double) {
            xd[0] = sigma * (x[1] - x[0]);
            xd[1] = R * x[0] - x[1] - x[0] * x[2];
            xd[2] = -b * x[2] + x[0] * x[1];
        }
    };

    /**
    * \brief a simple nonlinear first order differential equation
    **/
    struct ODE {
        void operator()(const std::vector<double>& x, std::vector<double>& xd, const double t) {
            xd[0] = cos(t);
        }
    };

    /**
    * \brief mass (https://en.wikipedia.org/wiki/Mass-spring-damper_model)
    **/
    struct Mass {
        Mass(std::vector<double>& state) : s(state), v(state) {}
    
        Simulation::StateParameter<double> s; // position
        Simulation::StateParameter<double> v; // velocity
        double m{}; // mass (zero mass means zero velocity)
        double f{}; // outside force
    
        void operator()(const std::vector<double>&, std::vector<double>& D, const double) {
            s(D) = v;
            v(D) = (m > 0.0) ? f / m : 0.0;
            f = 0.0;
        }
    };
    
    /**
    * \brief damper (https://en.wikipedia.org/wiki/Mass-spring-damper_model)
    **/
    struct Damper {
        Damper(Mass& b0, Mass& b1) : b0(b0), b1(b1) {}
    
        Mass& b0;
        Mass& b1;
    
        double dv{}; // velocity difference
        double c{}; // damping coefficient
        double f{}; // force
    
        void operator()(const std::vector<double>&, std::vector<double>&, const double) {
            dv = b0.v - b1.v;
            f = c * dv;
    
            b0.f -= f;
            b1.f += f;
        }
    };
    
    /**
    * \brief spring (https://en.wikipedia.org/wiki/Mass-spring-damper_model)
    **/
    struct Spring {
        Spring(Mass& b0, Mass& b1) : b0(b0), b1(b1) {
            l0 = b1.s - b0.s;
        }
    
        Mass& b0;
        Mass& b1;
    
        double l0{}; // initial spring length (distance between masses)
        double ds{}; // spring compression/extension
        double k{}; // spring coefficient
        double f{}; // force
    
        void operator()(const std::vector<double>&, std::vector<double>&, const double) {
            ds = l0 + b0.s - b1.s;
            f = k * ds;
    
            b0.f -= f;
            b1.f += f;
        }
    };
};

void test_airy() {
    Models::Airy system;
    std::vector<double> x{ { 1.0, 0.0 } };
    double t{};
    double dt{ 0.1 };
    double t_end{ 10.0 };

    Simulation::RK4 integrator;
    Simulation::History recorder;
    while (t < t_end) {
        recorder({ t, x[0], x[1] });
        integrator(system, x, t, dt);
    }

    recorder.csv("airy", { "t", "x0", "x1" }); // generate a file of comma separated values
}

void test_lorenz() {
    Models::Lorenz system;
    std::vector<double> x{ { 10.0, 1.0, 1.0 } };
    double t{};
    double dt{ 0.05 };
    double t_end{ 10.0 };

    Simulation::RK4 integrator;
    Simulation::History recorder;
    while (t < t_end) {
        recorder({ t, x[0], x[1], x[2] });
        integrator(system, x, t, dt);
    }

    recorder.csv("lorenz", { "t", "x0", "x1", "x2" });
};

void test_spring_damper() {
    std::vector<double> x;
    double t{};
    double dt{ 0.05 };
    double t_end{ 1.5 };

    // a spring and a damper are connecting one mass to another
    Models::Mass b0(x);
    Models::Mass b1(x);
    b1.m = 1.0;
    b1.s = 1.0;
    b1.v = 40.0;

    Models::Spring spring(b0, b1);
    spring.k = 2000.0;

    Models::Damper damper(b0, b1);
    damper.c = 5.0;

    Simulation::RK4 integrator;
    Simulation::History recorder;
    auto system = [&](const std::vector<double>& x, std::vector<double>& D, const double t) {
            spring(x, D, t);
            damper(x, D, t);
            b1(x, D, t);
    };

    while (t < t_end) {
        recorder({ t, b1.s });
        integrator(system, x, t, dt);
    }

    recorder.csv("spring-damper", { "t", "b1 position" });
};


void test_sampling() {
    std::vector<double> x{ { 0.0 } };
    double t{};
    double dt{ 0.1 };
    double t_end{ 10.0 };

    Models::ODE system;
    Simulation::RK4 integrator;
    Simulation::History recorder;
    while (t < t_end) {
        Simulation::SampleTime<double> sampler(t, dt);

        // evaluate every 0.25 seconds and trigger event when reaching 6.6 seconds
        if (sampler(0.25) || sampler.event(0.617)) {
            recorder({ t, x[0] });
        }

        integrator(system, x, t, dt);
    }

    recorder.csv("sampling", { "t", "x0" }); // generate a file of comma separated values
}

int main() {
    test_airy();
    test_lorenz();
    test_spring_damper();
    test_sampling();
}
