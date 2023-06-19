// lots of cpp11 code, no mpi code
#include <memory>
#include "dustmpi.h"
#include <dust/random/random.hpp>

using rng_state_type = dust::random::xoshiro256plus;

class model {
public:
  struct pars {
    double alpha;
    double beta;
    double gamma;
    double dt;
    int freq;
  };

  model(std::shared_ptr<const pars> pars) : pars_(pars) {
  }

  void update(const int time, const double * state,
              double * state_next, rng_state_type& rng_state) {
    double S = state[0];
    double I = state[1];
    double R = state[2];

    double N = S + I + R;

    double p_SI = 1 - std::exp(-(pars_->beta) * I / N);
    double p_IR = 1 - std::exp(-(pars_->gamma));
    double p_RS = 1 - std::exp(-(pars_->alpha));


    double n_SI = dust::random::binomial<double>(rng_state, S, p_SI * pars_->dt);
    double n_IR = dust::random::binomial<double>(rng_state, I, p_IR * pars_->dt);
    double n_RS = dust::random::binomial<double>(rng_state, R, p_RS * pars_->dt);

    state_next[0] = S - n_SI + n_RS;
    state_next[1] = I + n_SI - n_IR;
    state_next[2] = R + n_IR - n_RS;
    state_next[3] = (time % pars_->freq == 0) ? n_SI : state[3] + n_SI;
  }

private:
  std::shared_ptr<const pars> pars_;
};

std::vector<double> run_simulation(const model::pars& pars,
                                   const std::vector<double>& initial_state,
                                   const int end_time,
                                   const int n_particles,
                                   dust::random::prng<rng_state_type> * rng,
                                   bool use_mpi);
