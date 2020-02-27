#pragma once
#include <vector>
#include <fstream>
#include "array.h"
#include "velocity_field.h"
#include "level_set.h"
#include "visualizer.h"

template<int N>
class Simulator;

struct Material {
	typedef double real;
	real density; // negative for solids, 0 for free surface and positive for fluids.
	real surfaceTension;
	Material(real d, real s) : density(d), surfaceTension(s) {}
	~Material() {}
};

template<int N>
struct StaticCondition {
	// Environment
	typedef double real;
	StaticCondition() {};
	StaticCondition(const real* g, const std::vector<Material> &m) {
		copy_indexable<N>(g, gravity);
		materials = m;
	}

	real gravity[N];
	// unify single/multiphase sim.
	std::vector<Material> materials;
	
	void save(std::ofstream &os) {
		os << "staticcondition-begin" << std::endl;
		for (int i = 0; i < N; i++) {
			os << " " << gravity[i] << std::endl;
		}
		os << materials.size() << std::endl;
		for (int i = 0; i < materials.size(); i++) {
			os << materials[i].density << " " << materials[i].surfaceTension << std::endl;
		}
		os << "staticcondition-end" << std::endl;
	}

	void load(std::ifstream &is) {
		std::string token;
		is >> token;
		if (token != "staticcondition-begin") {
			std::cout << "failed to load static condition" << std::endl;
			exit(0);
		}
		for (int i = 0; i < N; i++) {
			is >> gravity[i];
		}
		int n_materials;
		is >> n_materials;
		if (n_materials <= 0) {
			std::cout << "no materials loaded" << std::endl;
		}
		for (int i = 0; i < n_materials; i++) {
			real rho, sigma;
			is >> rho >> sigma;
			materials.emplace_back(rho, sigma);
		}
		is >> token;
		if (token != "staticcondition-end") {
			std::cout << "failed to load static condition" << std::endl;
			exit(0);
		}
	}
};

template<int N>
struct DynamicCondition {
	typedef double real;

	DynamicCondition() {
		time = 0;
	}

	DynamicCondition(const real &_t, const std::vector<LevelSet<N>> &_l) :
		time(_t), levelSet(_l) {
	}

	real time;
	int n_steps = 0;
	std::vector<LevelSet<N>> levelSet;

	std::vector<int> blob_material;
	std::vector<Array<N, real>> a_var;
	std::vector<VelocityField<N>> v_var;
	mutable void *extra;

	//mutable std::vector<Array<N_DIM, int>> ghost_index_d;
	mutable Eigen::Matrix<real, Eigen::Dynamic, 1> p_G_d;

	void save(std::ofstream &os) const {
		os << "dynamiccondition-begin" << std::endl;
		os << time << std::endl;
		os << n_steps << std::endl;
		os << a_var.size() << std::endl;
		for (auto &a : a_var) {
			a.save(os);
		}
		os << v_var.size() << std::endl;
		for (auto &v : v_var) {
			v.save(os);
		}
		os << levelSet.size() << std::endl;
		for (auto &l : levelSet) {
			l.save(os);
		}
		os << blob_material.size() << std::endl;
		for (auto i : blob_material) {
			os << i << std::endl;
		}
		os << "dynamiccondition-end" << std::endl;
	}
	void load(std::ifstream &is) {
		std::string token;
		is >> token;
		if (token != "dynamiccondition-begin") {
			std::cout << "failed to load dynamiccondition" << std::endl;
		}
		is >> time >> n_steps;
		int n_a_var;
		is >> n_a_var;
		for (int i = 0; i < n_a_var; i++) {
			a_var.push_back(Array<N, real>(is));
		}
		int n_v_var;
		is >> n_v_var;
		for (int i = 0; i < n_v_var; i++) {
			v_var.push_back(VelocityField<N>(is));
		}
		int n_levelset;
		is >> n_levelset;
		for (int i = 0; i < n_levelset; i++) {
			LevelSet<N> l;
			l.load(is);
			levelSet.push_back(l);
		}
		int n_blob_material;
		is >> n_blob_material;
		for (int i = 0; i < n_blob_material; i++) {
			int mid;
			is >> mid;
			blob_material.push_back(mid);
		}
		is >> token;
		if (token != "dynamiccondition-end") {
			std::cout << "failed to load dynamic condition" << std::endl;
		}
	}
};

template<int N>
struct SimulatorSpec {
	typedef double real;

	int dim[N];
	int size;
	typename real domain[N][2];

	typename real dt;
	real dx;

	real getDt(const DynamicCondition<N>&) { return dt; }
	std::function<DynamicCondition<N>(const StaticCondition<N>&, const DynamicCondition<N>&, const SimulatorSpec<N>&)> stepper;
	std::function<void(const DynamicCondition<N>&)> visualizer;

	SimulatorSpec() {};
	SimulatorSpec (real dm[][2], const int *resolution, const real &_dt, std::string sim_name = "default") {
		size = 1;
		for (int i = 0; i < N; i++) {
			domain[i][0] = dm[i][0];
			domain[i][1] = dm[i][1];
			dim[i] = resolution[i];
			size *= dim[i];
		}
		dx = (domain[0][1] - domain[0][0]) / resolution[0];
		dt = _dt;
	}

	void save(std::ofstream &os) {
		os << "simulatorspec-begin" << std::endl;
		for (int i = 0; i < N; i++) {
			os << " " << dim[i];
		}
		os << std::endl;
		os << size << std::endl;
		for (int i = 0; i < N; i++) {
			os << domain[i][0] << " " << domain[i][1] << std::endl;
		}
		os << dt << std::endl;
		os << dx << std::endl;
		os << "simulatorspec-end" << std::endl;
	}
	void load(std::ifstream &is) {
		std::string token;
		is >> token;
		if (token != "simulatorspec-begin") {
			std::cout << "failed to load simulator spec: " << token << std::endl;
		}
		for (int i = 0; i < N; i++) {
			is >> dim[i];
		}
		is >> size;
		for (int i = 0; i < N; i++) {
			is >> domain[i][0] >> domain[i][1];
		}
		is >> dt >> dx;
		is >> token;
		if (token != "simulatorspec-end") {
			std::cout << "failed to load simulator spec: " << token << std::endl;
		}
	}
};

template<int N>
class Simulator {
public:
	typedef double real;
	Simulator() {};
	Simulator(const StaticCondition<N>& sc, const DynamicCondition<N> &dc, const SimulatorSpec<N> &ss) :
		staticCondition(sc), dynamicCondition(dc), simulatorSpec(ss) { 
	}

	
	StaticCondition<N> staticCondition;
	DynamicCondition<N> dynamicCondition;
	SimulatorSpec<N> simulatorSpec;

	void run(real endTime);

	void save(std::ofstream &os) {
		staticCondition.save(os);
		dynamicCondition.save(os);
		simulatorSpec.save(os);
	}
	void load(std::ifstream &is) {
		staticCondition.load(is);
		dynamicCondition.load(is);
		simulatorSpec.load(is);
	}
};

template<int N>
inline void Simulator<N>::run(real endTime) {
	while (dynamicCondition.time < endTime) {
		//std::ofstream os("result/ws-may21/" + std::to_string(dynamicCondition.n_steps) + ".txt");
		//save(os);
		//os.close();
		std::cout << "step starts: t=" << dynamicCondition.time << std::endl;
		dynamicCondition = simulatorSpec.stepper(staticCondition, dynamicCondition, simulatorSpec);
		std::cout << dynamicCondition.n_steps << "-th step finished, t=" << dynamicCondition.time << std::endl;
	}
}


