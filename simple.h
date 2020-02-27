#pragma once
#include "Array.h"
#include "Simulator.h"
#include "visualizer.h"

// vanilla step/project

namespace step_simple {
	typedef double real;
	typedef Eigen::Matrix<real, Eigen::Dynamic, 1> VectorX;

	// vanilla step
	template<int N>
	DynamicCondition<N> step(const StaticCondition<N>& sc, const DynamicCondition<N>& dc, const SimulatorSpec<N> &ss) {
		// condition check 
		bool check = false;
		check = check || (dc.levelSet.size() != 1);
		if (check) std::cout << "step function not applicable" << std::endl;

		VelocityField<N> vn = dc.v.extrapolateVelocity(dc.levelSet[0]);
		vn.setZeroBoundary();

		VelocityField<N> vm = vn;
		for (int i = 0; i < N; i++) {
			vm[i] += sc.gravity[i] * ss.dt;
		}

		vm = vn.advectVelocity(vm, ss.dt);
		LevelSet<N> advected_levelset = dc.levelSet[0];
		advected_levelset.sdf = vn.advect(dc.levelSet[0].sdf, ss.dt);
		vm.setZeroBoundary();

		advected_levelset.compute();

		DynamicCondition<N> dc_new(dc.time + ss.dt, vm, { advected_levelset }, dc.p);
		auto [p, v] = project(sc, dc_new, ss);
		dc_new.p = p; dc_new.v = v;
		return dc_new;
	}

	// vanilla project
	template<int N>
	std::tuple<Array<N, real>, VelocityField<N>>
		project(const StaticCondition<N>& sc, const DynamicCondition<N> &dc, const SimulatorSpec<N> &ss) {
		using namespace range_operator;
		Array<N, int> cellType = dc.levelSet[0].sdf.map<int>([](real sd) {return (sd <= 0) ? 0 : 1; });
		VectorX b(ss.size);
		std::vector<real> div = dc.v.divergence().flatten();
		for (int i = 0; i < ss.size; i++) {
			b[i] = -div[i];
		}

		Eigen::SparseMatrix<real> A(ss.size, ss.size);
		std::vector<Eigen::Triplet<real>> nnz;
		real scale = ss.dt / sc.materials[0].density / ss.dx / ss.dx;

		for (auto[ia, ct] : enumerate(cellType)) {
			int id = cellType._getLocation(ia);
			switch (ct) {
			case 1:
				nnz.emplace_back(id, id, 1.0);
				b[id] = 0.0;
				break;
			case 0:
				std::vector<int> neighbors = cellType.getNeighborsOf(ia);
				for (int nb : neighbors) {
					if (cellType[nb] == 0) {
						nnz.emplace_back(id, id, scale);
						nnz.emplace_back(id, nb, -scale);
					} else if (cellType[nb] == 1) {
						nnz.emplace_back(id, id, scale);
					}
				}
				break;
			}
		}

		A.setFromTriplets(nnz.begin(), nnz.end());
		Eigen::ConjugateGradient<Eigen::SparseMatrix<real>, Eigen::Upper, Eigen::IncompleteCholesky<real, Eigen::Upper>> cg;
		cg.compute(A);
		printf("begin solving\n");
		VectorX x = cg.solve(b);
		printf("iterations: %d, error: %f\n", cg.iterations(), cg.error());

		Array<N, real> p(ss.dim);
		p.setData(std::vector<real>(x.data(), x.data() + ss.size));

		using underscore_as_semicolon::_;
		VelocityField<N> v = dc.v;
		scale = ss.dt / sc.materials[0].density / ss.dx;
		for (int i = 0; i < N; i++) {
			auto dp_i = p[R<N>(i, { 1, _ })] - p[R<N>(i, { 0, -1 })];
			v[i][R<N>(i, { 1,-1 })] -= dp_i * scale;
		}
		return { p, v };
	}
}