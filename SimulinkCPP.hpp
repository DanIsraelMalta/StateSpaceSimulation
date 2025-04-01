#pragma once
#include <type_traits>
#include <concepts>
#include <vector>
#include <fstream>
#include <string>

/**
* collection of objects required for easy simulation of differential systems
**/
namespace Simulation {

    /**
    * internals
    **/
    namespace helpers {
        // traits to identify std::vector
        template<typename>   struct is_vector                 : std::false_type{};
        template<typename T> struct is_vector<std::vector<T>> : std::true_type {};
        template<typename T> constexpr bool is_vector_v = is_vector<T>::value;
    };

    /**
    * \brief an object which keeps timed samples while it is in scope
    **/
    template<typename T = double>
        requires(std::is_floating_point_v<T>)
    struct SampleTime {
        static constexpr T tol{ static_cast<T>(1e-8) };

        explicit constexpr SampleTime(T& _t, T& _dt) noexcept : t(_t), dt(_dt) {}

        SampleTime(const SampleTime& other) = delete;
        SampleTime& operator=(const SampleTime& other) = delete;
        SampleTime(SampleTime&& other) = delete;
        SampleTime& operator=(SampleTime&&) = delete;

        /**
        * \brief define sample rate for evaluation
        * @param {T,    in}  sample rate
        * @param {bool, out} true if time stamp is on sample rate
        **/
        bool operator()(const T sample_rate) noexcept {
            // housekeeping
            const std::size_t previous_samples_count{ static_cast<size_t>((this->t + tol) / sample_rate) };
            const T next_sample_time{ static_cast<T>(previous_samples_count + 1) * sample_rate };

            // delta time
            if (next_sample_time < this->t + this->dt - tol) {
                this->dt = next_sample_time - t;
            }

            // output
            return (this->t - next_sample_time + sample_rate < tol);
        }

        /**
        * \brief notify / trigger an event when given time stamp is reached
        * @param {T,    in} time stamp for event to be triggered.
        * @param {bool, in} true if time stamp reached
        **/
        constexpr bool event(const T event_time) noexcept {
            if (event_time < this->t + this->dt - tol &&
                event_time >= this->t + tol) {
                this->dt = event_time - this->t;
            }

            // output
            return (std::abs(event_time - this->t) < tol);
        }

        // internals
        private:
            T& t;
            T& dt;
    };

    /**
    * \brief parameter inside a larger state.
    *        states are vectors of parameters.
    **/
    template<typename T = double>
        requires(std::is_floating_point_v<T>)
    struct StateParameter {

        /**
        * \brief construct state parameter inside a given state and hold a reference to it
        **/
        constexpr explicit StateParameter(std::vector<T>& c, const T x0 = T{}) : index(c.size()), x(take_last(c, x0)) {}

        StateParameter(const StateParameter& other) = delete;
        StateParameter& operator=(const StateParameter& other) = delete;
        StateParameter(StateParameter&& other) = delete;
        StateParameter& operator=(StateParameter&&) = delete;

        /**
        * @param {vector<T>, in}  state holding parameter
        * @param {T&,        out} reference to parameter within state
        **/
        T& operator()(std::vector<T>& xd) const noexcept { return xd[index]; };

        /**
        * @param {T&, out} reference to parameter
        **/
        operator T&() const noexcept { return this->x; }

        /**
        * \brief assign new value to parameter
        * @param {T,  in}  new parameter value
        * @param {T&, out} reference to parameter
        **/
        T& operator=(const T y) noexcept { this->x = y; return this->x; }

        // internals
        private:
            const std::size_t index;
            T& x;

            /**
            * \brief workaround to allow references to be held inside objects
            * @param {vector<T>, in}  state to which parameter will be added
            * @param {T,         in}  parameter initial value
            * @param {T&,        out} reference to parameter within state
            **/
            T& take_last(std::vector<T>& c, const T x0) {
                c.emplace_back(x0);
                return c.back();
            }
    };
   
    /**
    * \brief an object which holds the history of a set of variables of identical type and can export it in CSV format
    **/
    template<class T = double>
    struct History {
 
        /**
        * \brief add row to history 
        * @param {vector<T>, in} row
        **/
        constexpr void operator()(std::vector<T> row) {
            this->history.emplace_back(row);
        }

        /**
        * \brief export history in CSV format
        * @param {string,         in}  path and name of CSV file
        * @param {vector<string>, in}  variables names associated with each column
        * @param {bool,           out} true if export was successful, false otherwise
        **/
        bool csv(const std::string& file_name, const std::vector<std::string>& names) {
            std::ofstream file;
            file.open(file_name + ".csv");
            if (!file) {
                return false;
            }

            const std::size_t num_names{ names.size() };
            for (std::size_t i{}; i < num_names; ++i) {
                file << names[i].c_str();

                if (i < num_names - 1) {
                    file << ",";
                }
            }
            if (num_names > 0) {
                file << '\n';
            }

            const std::size_t num_states{ this->history.front().size() };
            for (auto& step : this->history) {
                for (std::size_t i{}; i < num_states; ++i) {
                    file << step[i];
                    if (i < num_states - 1) {
                        file << ",";
                    }
                }
                file << '\n';
            }

            // output
            return true;
        }

        // internals
        private:
            std::vector<std::vector<T>> history; // history
    };


    /**
    * \brief fourth order Runge-Kutta integrator (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
    **/
    template<class STATE = std::vector<double>>
        requires(helpers::is_vector_v<STATE>)
    struct RK4 {
        using T = typename STATE::value_type;
        static_assert(std::is_floating_point_v<T>);

        /**
        * \brief perform one Runge-Kutta integration step in which both state and time are updated
        * @param {callable, in    } function depicting system to be integrated, its signature is void(const STATE&, STATE&, const T)
        * @param {STATE&,   in|out} system state (to be updated)
        * @param {T&,       in|out} system time (to be updated)
        * @param {T,        in    } time step
        **/
        template<class FUNC>
            requires(std::is_invocable_v<FUNC, const STATE&, STATE&, const T>)
        void operator()(FUNC&& func, STATE& x, T& t, const T dt) {
            // housekeeping
            const std::size_t len{ x.size() };
            const T t0{ t };
            const T dt_2{ static_cast<T>(0.5) * dt };
            const T dt_6{ static_cast<T>(1.0 / 6.0) * dt };
            if (this->xd.size() < len) {
                this->xd.resize(len);
                this->xd_temp.resize(len);
            }
            this->x0 = x;

            // 1st iteration
            func(x0, this->xd, t);
            for (std::size_t i{}; i < len; ++i) {
                x[i] = std::fma(dt_2, this->xd[i], this->x0[i]);
            }
            t += dt_2;

            // 2nd iteration
            func(x, this->xd_temp, t);
            for (std::size_t i{}; i < len; ++i) {
                this->xd[i] += 2 * this->xd_temp[i];
                x[i] = std::fma(dt_2, this->xd_temp[i], this->x0[i]);
            }

            // 3rd iteration
            func(x, this->xd_temp, t);
            for (std::size_t i{}; i < len; ++i) {
                this->xd[i] += 2 * this->xd_temp[i];
                x[i] = std::fma(dt, this->xd_temp[i], this->x0[i]);
            }
            t = t0 + dt;

            // 4th iteration
            func(x, this->xd_temp, t);
            for (std::size_t i{}; i < len; ++i) {
                x[i] = std::fma(dt_6, this->xd[i] + this->xd_temp[i], this->x0[i]);
            }
        }

        // internals
        private:
            STATE x0, xd, xd_temp;
    };    
}
