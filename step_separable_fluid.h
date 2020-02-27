#pragma once
// #include "symmetry_check.h"
#include "Simulator.h"
#include "visualizer.h"

namespace step_separable_fluid {
	typedef double real;
	typedef Eigen::Matrix<real, Eigen::Dynamic, 1> VectorX;


	template<int N>
	class LCP_Particle : public Particle<N> {
	public:
		static const int FREE = 0, SQUEEZED = 1, MERGED = 2;
		inline static const real T_SQUEEZED_THRESHOLD = 0.1;

		int status = FREE;
		real t_squeezed = 0;
	public:
		LCP_Particle() : Particle<N>() {}
		LCP_Particle(const LCP_Particle<N>& p2) : Particle<N>(p2) {
			status = p2.status;
			t_squeezed - p2.t_squeezed;
		}
		~LCP_Particle() {}
	};



	template<int N>
	DynamicCondition<N> step(const StaticCondition<N>& sc, const DynamicCondition<N>& dc, const SimulatorSpec<N> &ss) {
		bool check = true;
		check = check && (dc.a_var.size() == 1);
		check = check && (dc.blob_material.size() == dc.levelSet.size() - 1);
		check = check && (dc.levelSet.size() == dc.v_var.size());
		if (!check) {
			std::cout << "check failed" << std::endl;
			std::cout << dc.a_var.size() << ", " << dc.blob_material.size() << ", " << dc.levelSet.size() << std::endl;
		}



		// auto dc = split(_dc);
		const int SOLID = dc.levelSet.size() - 1;
		// extrapolate velocities
		std::vector<VelocityField<N>> vm = dc.v_var;
		for (int i = 0; i < dc.v_var.size() - 1; i++) {
			// vm[i] = dc.v_var[i].extrapolateVelocity(dc.levelSet[i]);
			vm[i] = dc.v_var[i].extrapolate_v(dc.levelSet[i].SDF());
		}

		auto[cell_type, __] = getCellType_m(dc.levelSet);

		std::vector<Array<N, real>> l_var(dc.levelSet.size());
		for (auto i = 0; i < l_var.size(); i++) {
			l_var[i] = dc.levelSet[i].SDF();
		}
		auto surface_info = get_surface_info(vm, l_var, cell_type);

		// apply gravity
		std::vector<VelocityField<N>> vn = vm;
		for (int i = 0; i < dc.v_var.size() - 1; i++) {
			for (int d = 0; d < N; d++) {
				vn[i][d] = vm[i][d] + sc.gravity[d] * ss.dt;
			}
		}

		// advect velocities 
		for (int i = 0; i < dc.v_var.size() - 1; i++) {
			vn[i] = vm[i].advectVelocity_cubic(vn[i], ss.dt);
		}



		// advect level sets
		std::vector<LevelSet<N>> ls = dc.levelSet;

		auto l_var1 = global_advect(l_var, surface_info, ss.dt, ss.dx);

		//std::cout << surface_info.size() << std::endl;
		//for (int i = 0; i < surface_info.size(); i++) {
		//	std::cout << i << std::endl;
		//	for (auto p : surface_info[i]) {
		//		std::cout << p[0] << ", " << p[1] << std::endl;
		//	}
		//}

		//int z = 10;
		//Visualizer vz(ss.dim[0] * z, ss.dim[1] * z, z);
		//vz.scalar(l_var[0]);
		//vz.mesh(l_var[0]);
		//vz.write("result/ws-may21/d" + std::to_string(dc.n_steps) + ".png");
		//vz.clear();
		//vz.scalar(l_var1[0]);
		//vz.mesh(l_var1[0]);
		//vz.write("result/ws-may21/e" + std::to_string(dc.n_steps) + ".png");
		//// exit(0); 

		for (int i = 0; i < dc.v_var.size() - 1; i++) {
			// ls[i].sdf = vm[i].advect_cubic(ls[i].sdf, ss.dt);
			ls[i].SDF() = l_var1[i];
			ls[i] = ls[i] - dc.levelSet[SOLID];
		}

		DynamicCondition<N> dc_new = dc;
		dc_new.v_var = vn;
		dc_new.levelSet = ls;

		// project
		auto dc_final = project_new(sc, dc_new, ss);
		dc_final.time = dc.time + ss.dt;
		dc_final.n_steps = dc.n_steps + 1;
		
		std::cout << dc_final.a_var.size() << dc_final.a_var[0].size() << std::endl;
		ss.visualizer(dc_final);
		// exit(0);
		
		if(0){
			auto &dc1 = dc_final;
			int l0 = 55, r0 = 95;
			int l1 = 25, r1 = 50;
			std::vector<std::vector<int>> r = { { l0, r0 },{ l1, r1 } };


			int z = 100;
			Visualizer vz((r0 - l0) * z, (r1 - l1) * z, z);

			Array<2, double> sp = dc1.a_var[0][r];

			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			vz.numbers(sp, { 0, 0 }, { -1.0 });
			vz.numbers(dc1.a_var[2][r], { 0, -50 }, { -1 });
			vz.numbers(dc1.a_var[1][r], { -50, 15 }, { -1 });
			vz.write("result/pressure.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			vz.numbers(dc1.a_var[5][r], { 0, 0 }, { 0 }, { 255, 255, 255 }, "%.0f");
			vz.numbers(dc1.a_var[7][r], { 0, -50 }, { 0 }, { 255, 255, 255 }, "%.0f");
			vz.numbers(dc1.a_var[6][r], { -50, 10 }, { 0 }, { 255, 255, 255 }, "%.0f");
			vz.write("result/id.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);
			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			auto &dc0 = dc;
			vz.numbers(dc0.v_var[1][0][r], { -50, -15 }, { 0.0 }, { 0, 255, 255 });
			vz.numbers(dc0.v_var[0][0][r], { -50, 15 }, { 0.0 }, { 255, 255, 0 });
			vz.numbers(dc0.v_var[1][1][r], { 0, -65 }, { 0.0 }, { 0, 255, 255 });
			vz.numbers(dc0.v_var[0][1][r], { 0, -35 }, { 0.0 }, { 255, 255, 0 });
			vz.write("result/v0.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);
			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			vz.numbers(dc1.v_var[1][0][r], { -50, -15 }, { 0.0 }, { 0, 255, 255 });
			vz.numbers(dc1.v_var[0][0][r], { -50, 15 }, { 0.0 }, { 255, 255, 0 });
			vz.numbers(dc1.v_var[1][1][r], { 0, -65 }, { 0.0 }, { 0, 255, 255 });
			vz.numbers(dc1.v_var[0][1][r], { 0, -35 }, { 0.0 }, { 255, 255, 0 });
			vz.write("result/velocity.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			Array<2, double> sl[3], cv[3];
			for (int i = 0; i < dc1.levelSet.size(); i++) {
				sl[i] = dc1.levelSet[i].SDF()[r];
				vz.mesh(sl[i]);
				cv[i] = getCurvature(sl[i]) / ss.dx;
			}

			vz.numbers(sl[1], { 0, -15 }, { 65535, -65535 }, { 0, 255, 255 });
			vz.numbers(sl[0], { 0, 15 }, { 65535, -65535 }, { 255, 255, 0 });
			vz.write("result/levelSet.png");

			vz.clear();
			vz.grid(r0 - l1, r1 - l1);
			Array<2, double> sl0[3], cv0[3];
			for (int i = 0; i < dc0.levelSet.size(); i++) {
				sl0[i] = dc0.levelSet[i].SDF()[r];
				vz.mesh(sl0[i]);
				cv[i] = getCurvature(sl0[i]) / ss.dx;
			}

			vz.numbers(sl0[1], { 0, -15 }, { 65535, -65535 }, { 0, 255, 255 });
			vz.numbers(sl0[0], { 0, 15 }, { 65535, -65535 }, { 255, 255, 0 });
			vz.write("result/l0.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			for (int i = 0; i < dc1.levelSet.size(); i++) {
				vz.mesh(sl[i]);
			}

			vz.numbers(cv[1], { 0, -15 }, { 0.0 }, { 0, 255, 255 });
			vz.numbers(cv[0], { 0, 15 }, { 0.0 }, { 255, 255, 0 });
			vz.write("result/curvature.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			double scale = 1;
			vz.numbers(Array<2, real>(dc1.a_var[4][r]), { 0, -50 }, { 65535 });
			vz.numbers(Array<2, real>(dc1.a_var[3][r]), { -50, 15 }, { 65535 });
			vz.write("result/gap.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			vz.numbers(Array<2, real>(dc0.a_var[4][r]), { 0, -50 }, { 65535 });
			vz.numbers(Array<2, real>(dc0.a_var[3][r]), { -50, 15 }, { 65535 });
			vz.write("result/gap0.png");

			vz.clear(); 
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			vz.points(surface_info[0], surface_info[3], { 0, 0 }, { l0, l1 });
			vz.points(surface_info[1], surface_info[4], { 0, 0 }, { l0, l1 }, { 0, -20 });
			vz.write("result/sp.png");

			vz.clear();
			vz.scalar(sp);
			vz.grid(r0 - l1, r1 - l1);

			for (auto &l : dc1.levelSet) {
				Array<2, double> sl = l.SDF()[r];
				vz.mesh(sl);
			}
			vz.indices(sp.dim, { l0, l1 });
			vz.write("result/indices.png");

			exit(0);
		}
		//if (dc.n_steps == 1)
		//	exit(0); 
		return dc_final;
	}

	// first half are surface points, second half are velocties
	template<int N>
	std::vector<std::vector<std::array<real, N>>> get_surface_info(const std::vector<VelocityField<N>> &v_var, std::vector<Array<N, real>> l_var, const Array<N, int> &cell_type) {
		using namespace range_operator;
		real TOL = 1e-6;
		std::vector<std::vector<std::array<real, N>>> ret;
		real band = 5.0;
		for (int i = 0; i < l_var.size(); i++) {
			auto sf = getSignFlip(l_var[i]);
			auto sp = initialize_b(l_var[i], sf);
			auto i2sp = gi2sp(l_var[i], sp, band);
			ret.push_back(sp);
		}
		for (int i = 0; i < l_var.size(); i++) {
			auto &sp = ret[i];
			auto vsp = sp;
			for (int j = 0; j < sp.size(); j++) {
				int ia[N];
				for (int k = 0; k < N; k++) {
					ia[k] = sp[j][k] + 0.5; // round to nearest cell center;
				}
				if (cell_type[ia] == i || !(0 <= cell_type[ia] && cell_type[ia] < l_var.size())) {
					vsp[j] = v_var[i].interpolate(sp[j]);
				} else {
					bool bounded = false;
					for (int k = 0; k < N; k++) {
						if (ia[k] > 0) {
							if (cell_type[I<N>(ia, k, -1)] == i) {
								if (abs(v_var[i][k][ia] - v_var[cell_type[ia]][k][ia]) < TOL ) {
									bounded = true;
								}
							}
						}
						if (ia[k] < l_var[0].dim[k] - 1) {
							if (cell_type[I<N>(ia, k, 1)] == i) {
								if (abs(v_var[i][k][I<N>(ia, k, 1)] - v_var[cell_type[ia]][k][I<N>(ia, k, 1)]) < TOL) {
									bounded = true;
								}
							}
						}
					}
					if (bounded) {
						vsp[j] = v_var[cell_type[ia]].interpolate(sp[j]);
					} else {
						vsp[j] = v_var[i].interpolate(sp[j]);
					}
				}
			}
			ret.push_back(vsp);
		}
		return ret;
	}

	template<int N>
	std::vector<Array<N, real>> global_advect(std::vector<Array<N, real>> l_var, const std::vector<std::vector<std::array<real, N>>> &surface_info, real dt, real dx) {
		real band = 5.0;
		std::vector<Array<N, int>> i2sp(l_var.size());
		for (int i = 0; i < l_var.size(); i++) {
			auto sf = getSignFlip(l_var[i]);
			auto sp = initialize_b(l_var[i], sf);
			i2sp[i] = gi2sp(l_var[i], sp, band);
		}

		auto ret = l_var;
		for (int i = 0; i < l_var.size(); i++) {
			ret[i] = advect(l_var[i], surface_info[i], surface_info[i + l_var.size()], i2sp[i], dt, dx);
		}
		return ret;
	}

	template<int N>
	std::vector<VelocityField<N>> global_extrapolation(const std::vector<VelocityField<N>> &v_var, const Array<N, int> &cell_type) {
		using namespace range_operator;

		auto ret = v_var;

		const int AIR = v_var.size();
		const int VALID = 1;
		std::vector<Array<N, int>> tag(v_var.size(), Array<N, int>(cell_type.dim, 0));
		for (auto[ia, ct] : enumerate(cell_type)) {
			if (ct != AIR) {
				tag[ct][ia] = VALID;
			}
		}



		std::vector<std::array<Array<N, int>, N>> valid(AIR);
		for (int ct = 0; ct < AIR; ct++) {
			for (int i = 0; i < N; i++) {
				valid[ct][i] = Array<N, int>(v_var[ct][i].dim, 0);
			}
		}
		for (auto[ia, ct] : enumerate(cell_type)) {
			if (ct == AIR) {
				continue;
			}
			for (int i = 0; i < N; i++) {
				valid[ct][i][ia] = 1;
				valid[ct][i][I<N>(ia, i, 1)] = 1;
			}
		}
/*
		int z = 20;
		Visualizer vz(cell_type.dim[0] * z, cell_type.dim[1] * z, z);
		vz.clear();
		// vz.scalar(tag[0]);
		vz.grid(cell_type.dim[0], cell_type.dim[1]);
		vz.scalar_dots(valid[0][0], { -z / 2, 0 });
		vz.scalar_dots(valid[0][1], { 0, -z / 2 });
		vz.write("result/0.png");
		vz.clear();
		// vz.scalar(tag[1]);
		vz.grid(cell_type.dim[0], cell_type.dim[1]);
		vz.scalar_dots(valid[1][0], { -z / 2, 0 });
		vz.scalar_dots(valid[1][1], { 0, -z / 2 });
		vz.write("result/1.png");
		vz.clear();
		// vz.scalar(tag[2]);
		vz.grid(cell_type.dim[0], cell_type.dim[1]);
		vz.scalar_dots(valid[2][0], { -z / 2, 0 });
		vz.scalar_dots(valid[2][1], { 0, -z / 2 });
		vz.write("result/2.png");
		*/
		auto eq = [](real a, real b) { return abs(a - b) < 1e-4; };

		for (int ct = 0; ct < AIR; ct++) {
			for (auto ia : Range<N>(tag[ct].dim)) {
				if (tag[ct][ia] == VALID) {
					for (int i = 0; i < N; i++) {
						if (ia[i] > 0) {
							int ct_l = cell_type[I<N>(ia, i, -1)];
							if (ct_l != ct && ct_l != AIR) {
								if (eq(ret[ct_l][i][ia], ret[ct][i][ia])) {
									auto e_ia = I<N>(ia, i, -1);
									for (int j = 0; j < N; j++) {
										auto e_ia_j_r = I<N>(e_ia, j, 1);
										if (valid[ct][j][e_ia] != 1) {
											ret[ct][j][e_ia] = ret[ct_l][j][e_ia];
											valid[ct][j][e_ia] = 2;
										}
										if (valid[ct][j][e_ia_j_r] != 1) {
											ret[ct][j][e_ia_j_r] = ret[ct_l][j][e_ia_j_r];
											valid[ct][j][e_ia_j_r] = 2;
										}
									}
								}
							}
						}
						if (ia[i] + 1 < tag[ct].dim[i]) {
							int ct_r = cell_type[I<N>(ia, i, 1)];
							if (ct_r != ct && ct_r != AIR) {
								if (eq(ret[ct_r][i][I<N>(ia, i, 1)], ret[ct][i][I<N>(ia, i, 1)])) {
									auto e_ia = I<N>(ia, i, 1);
									for (int j = 0; j < N; j++) {
										auto e_ia_j_r = I<N>(e_ia, j, 1);
										if (valid[ct][j][e_ia] != 1) {
											ret[ct][j][e_ia] = ret[ct_r][j][e_ia];
											valid[ct][j][e_ia] = 2;
										}
										if (valid[ct][j][e_ia_j_r] != 1) {
											ret[ct][j][e_ia_j_r] = ret[ct_r][j][e_ia_j_r];
											valid[ct][j][e_ia_j_r] = 2;
										}
									}
								}
							}
						}
					}
				}
			}
		}
/*
		vz.grid(cell_type.dim[0], cell_type.dim[1]);
		vz.scalar_dots(valid[0][0], { -z / 2, 0 });
		vz.scalar_dots(valid[0][1], { 0, -z / 2 });
		vz.write("result/00.png");
		vz.clear();
		// vz.scalar(tag[1]);
		vz.grid(cell_type.dim[0], cell_type.dim[1]);
		vz.scalar_dots(valid[1][0], { -z / 2, 0 });
		vz.scalar_dots(valid[1][1], { 0, -z / 2 });
		vz.write("result/10.png");
		vz.clear();
		// vz.scalar(tag[2]);
		vz.grid(cell_type.dim[0], cell_type.dim[1]);
		vz.scalar_dots(valid[2][0], { -z / 2, 0 });
		vz.scalar_dots(valid[2][1], { 0, -z / 2 });
		vz.write("result/20.png");
*/
		for (int ct = 0; ct < AIR; ct++) {
			for (int i = 0; i < N; i++) {
				extrapolate(valid[ct][i], ret[ct].v[i], 30);
			}
		}

		return ret;
	}


	// return false if there is no conflict
	int resolve_conflict(real *l, int n) {
		real MAX = 65535;
		real min[2] = { MAX, MAX };
		int t[2] = { -1, -1 };
		for (int i = 0; i < n; i++) {
			if (l[i] <= 1) {
				if (l[i] <= min[0]) {
					min[1] = min[0];
					t[1] = t[0];
					min[0] = l[i];
					t[0] = i;
				} else if (l[i] <= min[1]) {
					min[1] = l[i];
					t[1] = i;
				}
			}
		}
		if (min[0] > 0) { // smallest sdf > 0: outside every levelset
			return n + 1;
		} else {
			real mean = (min[0] + min[1]) / 2.0;
			if (mean <= 0) {
				for (int i = 0; i < n; i++) {
					l[i] -= mean;
				}
			}
			if (min[0] == min[1]) {
				l[t[0]] -= 1e-13;
			}
			return t[0];
		}
	}

	int tweak(real *l, int n) {
		real min[2] = { 1, 1 };
		int t[2] = { -1, -1 };
		for (int i = 0; i < n; i++) {
			if (l[i] <= 0.5) {
				if (l[i] <= min[0]) {
					min[1] = min[0];
					t[1] = t[0];
					min[0] = l[i];
					t[0] = i;
				} else if (l[i] <= min[1]) {
					min[1] = l[i];
					t[1] = i;
				}
			}
		}
		if (t[0] < 0) {
			return n + 1;
		} else if (t[1] < 0) {
			return n + 1;
		} else {
			if (min[0] + min[1] < 0.5) {
				l[t[0]] -= min[0] + 1e-14;
				l[t[0]] += min[0] + 1e-14;
			}
			return t[0];
		}
	}

	// when more than 2 fluid levelsets are present
	template<int N>
	std::tuple<Array<N, int>, std::vector<LevelSet<N>>> getCellType_m(std::vector<LevelSet<N>> levelSet) {
		using namespace range_operator;
		Array<N, int> cellType(levelSet[0].SDF().dim);
		int n_fluids = levelSet.size() - 1;
		real * l = new real[n_fluids];
		for (auto &[ia, ct] : enumerate(cellType)) {
			if (levelSet[n_fluids][ia] < 0) {
				ct = n_fluids; // solid
				continue;
			}
			for (int i = 0; i < n_fluids; i++) {
				l[i] = levelSet[i][ia];
			}
			ct = resolve_conflict(l, n_fluids);
			for (int i = 0; i < n_fluids; i++) {
				levelSet[i][ia] = l[i];
			}

		}
		for (int i = 0; i < n_fluids; i++) {
			levelSet[i].redistance();
		}
		return { cellType, levelSet };
	}


	// to be modified to m version
	template<int N>
	std::array<std::tuple<int, real>, 2 * N + 2> getVelocityUpdateFormula_m(int fluid_type, const int *ia, int axis, const Array<N, int> &cellType, const Array<N, int> &id,
		const std::vector<LevelSet<N>> &levelset, const Array<N, int>* ghostId, const std::vector<VelocityField<N>> &v_var, const std::vector<Array<N, real>> &jump,
		const StaticCondition<N> &sc, const SimulatorSpec<N> &ss, const std::vector<int> &blob_material) {

		typedef std::tuple<int, real> pair;

		int blob_id = fluid_type;
		int n_fluids = blob_material.size();
		const int SOLID = n_fluids;
		const int AIR = SOLID + 1;

		if (blob_id == SOLID) // solid
			return { pair{ 0, 0 }, pair{ -1, 0.0 } };

		using namespace range_operator;

		for (int i = 0; i < N; i++) {
			if (ia[i] <= 0 || ia[i] >= id.dim[i]) {
				return { pair{ 0, 0.0 }, pair{ -1, 0.0 } };
			}
		}

		auto isFluid = [n_fluids](int t) {return t < n_fluids; };

		auto ia_l = I<N>(ia, axis, -1);
		auto ia_r = ia;
		int ct_l = cellType[ia_l];
		int ct_r = cellType[ia_r];

		int blob_id_l = ct_l;
		int blob_id_r = ct_r;

		std::array<pair, 2 * N + 2> zero = { pair{ 0, 0.0 }, pair{ -1, 0.0 } };

		int material_id = blob_material[blob_id];

		real scale = ss.dt / sc.materials[material_id].density / ss.dx;
		real u_star = v_var[blob_id][axis][ia];
		real sdf_l = levelset[blob_id].SDF()[ia_l];
		real sdf_r = levelset[blob_id].SDF()[ia_r];

		if (blob_id_l == blob_id) {
			real J = jump[blob_id][ia_l];
			int id_l = id[ia_l];
			real theta = sdf_l / (sdf_l - sdf_r);
			theta = std::clamp<real>(theta, 1e-3, 1 - 1e-3);

			if (blob_id_r == blob_id) {
				int id_r = id[ia_r];
				return { pair{ 0, u_star }, pair{ id_l, scale }, pair{ id_r, -scale }, pair{ -1, 0.0 } };
			} else if (isFluid(blob_id_r)) {
				int id_r = ghostId[axis][ia];
				return { pair{ 0, u_star - J * scale / theta }, pair{ id_l, scale / theta }, pair{ id_r, -scale / theta }, pair{ -1, 0.0 } };
			} else if (blob_id_r == AIR) {
				return { pair{ 0, u_star - J * scale / theta }, pair{ id_l, scale / theta }, pair{ -1, 0.0 } };
			} else if (blob_id_r == SOLID) {
				int id_r = ghostId[axis][ia];
				return { pair{ 0, u_star - J * scale / theta }, pair{ id_l, scale / theta }, pair{ id_r, -scale / theta }, pair{ -1, 0.0 } };
			}
		} else if (blob_id_r == blob_id) {
			real J = jump[blob_id][ia_r];
			int id_r = id[ia_r];
			real theta = sdf_r / (sdf_r - sdf_l);
			theta = std::clamp<real>(theta, 1e-3, 1 - 1e-3);

			if (isFluid(blob_id_l)) {
				// different fluid on the left
				int id_l = ghostId[axis][ia];
				return { pair{ 0, u_star + J * scale / theta }, pair{ id_l, scale / theta }, pair{ id_r, -scale / theta }, pair{ -1, 0.0 } };
			} else if (blob_id_l == AIR) {
				return { pair{ 0, u_star + J * scale / theta }, pair{ id_r, -scale / theta }, pair{ -1, 0.0 } };
			} else if (blob_id_l == SOLID) {
				int id_l = ghostId[axis][ia];
				return { pair{ 0, u_star + J * scale / theta }, pair{ id_l, scale / theta }, pair{ id_r, -scale / theta }, pair{ -1, 0.0 } };
			}
		}

		return zero;
	}

	// output: id, ghostId: n_cell, n_ghost_cell
	template<int N>
	std::tuple<int, int> getId(Array<N, int> cellType, Array<N, int>& id, Array<N, int> *ghostId, const std::vector<LevelSet<N>>& levelSet, std::vector<real> &ghostGap) {
		const int n_fluids = levelSet.size() - 1;
		auto isFluid = [n_fluids](int t) {return t < n_fluids; };
		const int SOLID = n_fluids;
		using namespace range_operator;
		int cnt = 1;
		for (auto[ia, ct] : enumerate(cellType)) {
			if (isFluid(ct)) {
				id[ia] = cnt;
				cnt++;
			}
		}
		int cnt_G = cnt;
		VelocityField<N> v(cellType.dim, 1);
		ghostGap = std::vector<real>(cnt, 65535);
		for (int i = 0; i < N; i++) {
			for (auto ia : Subrange<N>(v[i].dim, R<N>(i, { 1, -1 }))) {
				int ct_l = cellType[I<N>(ia, i, -1)];
				int ct_r = cellType[ia];
				bool fluid_fluid = isFluid(ct_l) && isFluid(ct_r) && ct_l != ct_r;
				bool fluid_solid = (isFluid(ct_l) && ct_r == SOLID) || (isFluid(ct_r) && ct_l == SOLID);
				if (fluid_fluid || fluid_solid) {
					ghostId[i][ia] = cnt_G;
					cnt_G++;

					auto ia_l = I<N>(ia, i, -1);
					auto ia_r = ia;
					int ct_l = cellType[I<N>(ia, i, -1)];
					int ct_r = cellType[ia];
					real sdf0_l = levelSet[ct_l].SDF()[ia_l];
					real sdf0_r = levelSet[ct_l].SDF()[ia_r];
					real sdf1_l = levelSet[ct_r].SDF()[ia_l];
					real sdf1_r = levelSet[ct_r].SDF()[ia_r];
					real theta0 = sdf0_l / (sdf0_l - sdf0_r);
					real theta1 = sdf1_r / (sdf1_r - sdf1_l);
					theta0 = std::clamp<real>(theta0, 1e-3, 1 - 1e-3);
					theta1 = std::clamp<real>(theta1, 1e-3, 1 - 1e-3);

					ghostGap.push_back(std::max<real>(0, 1 - theta0 - theta1));
				}
			}
		}
		return { cnt, cnt_G - cnt };
	}

	template<int N>
	std::tuple<Eigen::SparseMatrix<real>, VectorX> getMatrix_m(int sz, const Array<N, int> &cellType, const Array<N, int> &id, const Array<N, int> *ghostId,
		const std::vector<std::array<Array<N + 1, int>, N>> &uf_id, const std::vector<std::array<Array<N + 1, real>, N>> &uf_coeff,
		const std::vector<real> &ghostGap, const SimulatorSpec<N> &ss) {
		using namespace range_operator;

		const int SOLID = uf_id.size() - 1;
		const int AIR = SOLID + 1;
		auto isFluid = [SOLID](int t) {return t < SOLID; };

		VectorX b(sz);
		b = Eigen::Matrix<real, Eigen::Dynamic, 1>::Zero(sz);

		std::vector<Eigen::Triplet<real>> nnz;
		nnz.emplace_back(0, 0, 1.0);
		b[0] = -1.0;

		for (auto &[ia, ct] : enumerate(cellType)) {
			// build equation for all fluid cells
			if (isFluid(ct)) {
				int self = id[ia];

				for (int i = 0; i < N; i++) {

					for (int k = 0; k < 2 * N + 2; k++) {
						int _id = uf_id[ct][i][A<N>(ia, k)];
						real _coeff = uf_coeff[ct][i][A<N>(ia, k)];
						if (_id < 0) break;
						if (_id == 0) {
							b[self] -= _coeff;
						} else {
							nnz.emplace_back(self, _id, -_coeff);
						}
					}

					for (int k = 0; k < 2 * N + 2; k++) {
						int _id = uf_id[ct][i][A<N>(I<N>(ia, i, 1), k)];
						real _coeff = uf_coeff[ct][i][A<N>(I<N>(ia, i, 1), k)];
						if (_id < 0) break;
						if (_id == 0) {
							b[self] += _coeff;
						} else {
							nnz.emplace_back(self, _id, _coeff);
						}
					}
				}
			}
		}

		for (int i = 0; i < N; i++) {
			for (auto[ia, g_id] : enumerate(ghostId[i])) {
				if (g_id > 0) {
					auto ia_l = I<N>(ia, i, -1);
					auto ia_r = ia;
					int ct_l = cellType[ia_l];
					int ct_r = cellType[ia_r];

					for (int k = 0; k < 2 * N + 2; k++) {
						int _id = uf_id[ct_l][i][A<N>(ia, k)];
						real _coeff = uf_coeff[ct_l][i][A<N>(ia, k)];
						if (_id < 0) break;
						if (_id == 0) {
							b[g_id] -= _coeff;
						} else {
							nnz.emplace_back(g_id, _id, -_coeff);
						}
					}


					for (int k = 0; k < 2 * N + 2; k++) {
						int _id = uf_id[ct_r][i][A<N>(ia, k)];
						real _coeff = uf_coeff[ct_r][i][A<N>(ia, k)];
						if (_id < 0) break;
						if (_id == 0) {
							b[g_id] += _coeff;
						} else {
							nnz.emplace_back(g_id, _id, _coeff);
						}
					}

					b[g_id] += ghostGap[g_id] * ss.dx / ss.dt;
				}
			}
		}
		Eigen::SparseMatrix<real> A(sz, sz);

		A.setFromTriplets(nnz.begin(), nnz.end());

		return { A, b };
	}


	template<int N>
	VelocityField<N> getVelocity(const SimulatorSpec<N> &ss, const std::array<Array<N + 1, int>, N> &uf_id, const std::array<Array<N + 1, real>, N> &uf_coeff, const VectorX &p) {
		using namespace range_operator;
		VelocityField<N> ret(ss.dim, ss.dx);
		for (int i = 0; i < N; i++) {
			ret[i] = 0;
			for (auto &[ia, id] : enumerate(uf_id[i])) {
				if (id < 0) continue;
				ret[i][ia] += p[id] * uf_coeff[i][ia];
			}
		}
		return ret;
	}

	template<int N>
	Array<N, real> getPressure(const SimulatorSpec<N> &ss, const Array<N, int> &id, const VectorX &p) {
		using namespace range_operator;
		Array<N, real> ret(ss.dim);
		ret = 0;
		for (auto &[ia, val] : enumerate(ret)) {
			if (id[ia] > 0) {
				val = p[id[ia]];
			}
		}
		return ret;
	}

	template<int N>
	std::tuple<Array<N, int>, std::vector<LevelSet<N>>> tweakLevelSet_m(Array<N, int> cellType, std::vector<LevelSet<N>> levelset) {
		// return { cellType, levelset };
		using namespace range_operator;
		const int n_fluids = levelset.size() - 1;
		const int AIR = n_fluids + 1;
		real * l = new real[n_fluids];

		for (auto &[ia, ct] : enumerate(cellType)) {
			if (ct != AIR)
				continue;
			for (int i = 0; i < n_fluids; i++) {
				l[i] = levelset[i][ia];
			}
			ct = tweak(l, n_fluids);
			for (int i = 0; i < n_fluids; i++) {
				levelset[i][ia] = l[i];
			}
		}

		return { cellType, levelset };
	}

	template<int N>
	std::vector<LCP_Particle<N>> reseedParticles(const Array<N, int> *ghostId, const VectorX &p, std::vector<LCP_Particle<N>> particles, const SimulatorSpec<N> &ss) {
		using namespace range_operator;
		using underscore_as_semicolon::_;
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<typename Particle<N>::real> dis(0.0, 1.0);

		std::vector<LCP_Particle<N>> ret;
		Array<N, int> si(ss.dim);
		si = 0;
		for (int i = 0; i < N; i++) {
			si = si + (ghostId[i][R<N>(i, { 0, -1 })] + ghostId[i][R<N>(i, { 1, _ })]);
		}
		Array<N, std::vector<int>> sorted = sortParticles(ss.dim, particles);
		for (auto ia : Range<N>(ss.dim)) {
			// if it is an LCP boundary, check whether neighboring cells have enough particles
			if (si[ia] > 0) {
				auto samples = sampleParticle<N, LCP_Particle>(ia, 0.5, 16, dis, gen);
				int status = LCP_Particle<N>::FREE;
				for (int i = 0; i < N; i++) {
					if (ghostId[i][ia] > 0) {
						if (p[ghostId[i][ia]] > 1e-4) {
							status = LCP_Particle<N>::SQUEEZED;
						}
					}
					if (ghostId[i][I<N>(ia, i, 1)] > 0) {
						if (p[ghostId[i][I<N>(ia, i, 1)]] > 1e-4) {
							status = LCP_Particle<N>::SQUEEZED;
						}
					}
				}
				if (status == LCP_Particle<N>::SQUEEZED) {
					real time = 0.0;
					for (auto i : sorted[ia]) {
						time = std::max<real>(time, particles[i].t_squeezed);
					}
					for (auto &s : samples) {
						s.t_squeezed = time;
						s.status = status;
					}
				}
				ret.insert(ret.end(), samples.begin(), samples.end());
			}
		}
		return ret;
	}

	template<int N>
	DynamicCondition<N> split(const DynamicCondition<N>& dc) {
		int n_fluids = dc.levelSet.size() - 1;
		DynamicCondition<N> ret = dc;
		ret.levelSet.clear();
		ret.v_var.clear();
		ret.blob_material.clear();
		for (int i = 0; i < n_fluids; i++) {
			if (!dc.levelSet[i].empty()) {
				int nc = dc.levelSet[i].identifyComponents();
				auto ln = dc.levelSet[i].split(nc);
				ret.levelSet.insert(ret.levelSet.end(), ln.begin(), ln.end());
				ret.v_var.insert(ret.v_var.end(), nc, dc.v_var[i]);
				ret.blob_material.insert(ret.blob_material.end(), nc, dc.blob_material[i]);
			}
		}
		ret.levelSet.push_back(dc.levelSet[n_fluids]);
		ret.v_var.push_back(dc.v_var[n_fluids]);
		ret.extra = dc.extra;
		return ret;
	}

	template<int N>
	DynamicCondition<N>
		project_new(const StaticCondition<N>& sc, const DynamicCondition<N> &dc, const SimulatorSpec<N> &ss) {

		using underscore_as_semicolon::_;
		using namespace range_operator;




		const int n_fluids = dc.levelSet.size() - 1;

		// setup cellType and ID
 		auto[cellType0, levelSet0] = getCellType_m(dc.levelSet);
		auto[cellType, levelSet] = tweakLevelSet_m(cellType0, levelSet0);


		Array<N, int> id(ss.dim), ghostId[N];
		id = 0;
		for (int i = 0; i < N; i++) {
			ghostId[i] = Array<N, int>(dc.v_var[0][i].dim);
			ghostId[i] = 0;
		}

		std::vector<real> ghostGap;
		auto[nCell, nGhost] = getId<N>(cellType, id, ghostId, levelSet, ghostGap);

		std::vector<std::array<Array<N + 1, int>, N>> uf_id(n_fluids + 1);
		std::vector<std::array<Array<N + 1, real>, N>> uf_coeff(n_fluids + 1);
		std::vector<Array<N, real>> jump(n_fluids);
		
		for (int i = 0; i < n_fluids; i++) {
			jump[i] = levelSet[i].curvature(ss.dx) * sc.materials[dc.blob_material[i]].surfaceTension;
			// jump[i] = 0;
			//real dt2 = ss.dt / 100.0;
			//auto temp = levelSet[i].stepMeanCurvatureMotion(ss.dx, sc.materials[dc.blob_material[i]].surfaceTension, dt2, 100);
			//jump[i] = (temp.sdf - levelSet[i].sdf) * (ss.dx / ss.dt);
			//jump[i] = jump[i] * getSignFlip(levelSet[i].sdf);
		}

		//for (int i = 0; i < n_fluids; i++) {
		//	jump[i] = Array<N_DIM, real>(ss.dim, 0.0);
		//}

		for (int i = 0; i < N; i++) {
			for (int id = 0; id < n_fluids + 1; id++) {
				uf_id[id][i] = Array<N + 1, int>(A<N>(dc.v_var[id][i].dim, 2 * N + 2));
				uf_id[id][i] = -1;
				uf_coeff[id][i] = Array<N + 1, real>(A<N>(dc.v_var[id][i].dim, 2 * N + 2));
			}
			for (int id = 0; id < n_fluids + 1; id++) {
				for (auto ia : Range<N>(dc.v_var[id][i].dim)) {
					auto coeff_pairs = getVelocityUpdateFormula_m<N>(id, ia, i, cellType, id, levelSet, ghostId, dc.v_var, jump, sc, ss, dc.blob_material);
					int k = 0;
					for (auto[_id, coeff] : coeff_pairs) {
						uf_id[id][i][A<N>(ia, k)] = _id;
						uf_coeff[id][i][A<N>(ia, k)] = coeff;
						k++;
						if (_id == -1) break;
					}
				}
			}
		}

		auto[A, b] = getMatrix_m<N>(nCell + nGhost, cellType, id, ghostId, uf_id, uf_coeff, ghostGap, ss);


		VectorX p = MLCP_Eigen::mlcp_solve(A, b, nCell, nGhost, 100);

		auto dc_new = dc;
		for (int i = 0; i < n_fluids; i++) {
			dc_new.v_var[i] = getVelocity<N>(ss, uf_id[i], uf_coeff[i], p);
		}

		// dump ghost values
		Array<N, real> ghost_values[N];
		Array<N, real> ghost_gaps[N];
		for (int i = 0; i < N; i++) {
			ghost_values[i] = Array<N, real>(ghostId[i].dim, -1);
			ghost_gaps[i] = Array<N, real>(ghostId[i].dim, 65535);
			for (int j = 0; j < ghostId[i].size(); j++) {
				if (ghostId[i][j] != 0) {
					ghost_values[i][j] = p[ghostId[i][j]];
					ghost_gaps[i][j] = ghostGap[ghostId[i][j]];
				}
			}
		}

		//auto flipped_id = sym_check::flip(id, 0);
		//auto flipped_gid0 = sym_check::flip(ghostId[0], 0);
		//auto flipped_gid1 = sym_check::flip(ghostId[1], 0);
		//VectorX flipped_p(p.size());
		//for (auto ia : Range<N_DIM>(ss.dim)) {
		//	if (flipped_id[ia] > 0) {
		//		if (!id[ia] > 0) {
		//			std::cout << "error" << std::endl;
		//		}
		//		if (abs(b[id[ia]] - b[flipped_id[ia]]) > 1e-8) {
		//			std::cout << id[ia] << ", " << flipped_id[ia] << ": " << b[id[ia]] << ", " << b[flipped_id[ia]] << std::endl;
		//		}
		//		flipped_p[flipped_id[ia]] = p[id[ia]];
		//	}
		//}

		//for (auto ia : Range<N_DIM>(ghostId[0].dim)) {
		//	if (flipped_gid0[ia] > 0) {
		//		if (!ghostId[0][ia] > 0) {
		//			std::cout << "error" << std::endl;
		//		}
		//		flipped_p[flipped_gid0[ia]] = p[ghostId[0][ia]];
		//	}
		//}

		//for (auto ia : Range<N_DIM>(ghostId[1].dim)) {
		//	if (flipped_gid1[ia] > 0) {
		//		if (!ghostId[1][ia] > 0) {
		//			std::cout << "error" << std::endl;
		//		}
		//		flipped_p[flipped_gid1[ia]] = p[ghostId[1][ia]];
		//	}
		//}
		//int ct = 0;
		//int ia[2] = { 76, 34 };
		//std::cout << uf_id[ct][0][{ia[0], ia[1], 0}] << ", " << uf_id[ct][0][{ia[0], ia[1], 1}] << ", " << uf_id[ct][0][{ia[0], ia[1], 2}] << ", " << uf_id[ct][0][{ia[0], ia[1], 3}] << "|" << std::endl;
		//std::cout << uf_coeff[ct][0][{ia[0], ia[1], 0}] << ", " << uf_coeff[ct][0][{ia[0], ia[1], 1}] << ", " << uf_coeff[ct][0][{ia[0], ia[1], 2}] << ", " << uf_coeff[ct][0][{ia[0], ia[1], 3}] << "|" << std::endl;
		//std::cout << uf_id[ct][1][{ia[0], ia[1], 0}] << ", " << uf_id[ct][1][{ia[0], ia[1], 1}] << ", " << uf_id[ct][1][{ia[0], ia[1], 2}] << "|" << std::endl;
		//std::cout << uf_coeff[ct][1][{ia[0], ia[1], 0}] << ", " << uf_coeff[ct][1][{ia[0], ia[1], 1}] << ", " << uf_coeff[ct][1][{ia[0], ia[1], 2}] << "|" << std::endl;

		//std::cout << uf_id[ct][0][{ia[0] + 1, ia[1], 0}] << ", " << uf_id[ct][0][{ia[0] + 1, ia[1], 1}] << ", " << uf_id[ct][0][{ia[0] + 1, ia[1], 2}] << ", " << uf_id[ct][0][{ia[0] + 1, ia[1], 3}] << "|" << std::endl;
		//std::cout << uf_coeff[ct][0][{ia[0] + 1, ia[1], 0}] << ", " << uf_coeff[ct][0][{ia[0] + 1, ia[1], 1}] << ", " << uf_coeff[ct][0][{ia[0] + 1, ia[1], 2}] << ", " << uf_coeff[ct][0][{ia[0] + 1, ia[1], 3}] << "|" << std::endl;
		//std::cout << uf_id[ct][1][{ia[0], ia[1] + 1, 0}] << ", " << uf_id[ct][1][{ia[0], ia[1] + 1, 1}] << ", " << uf_id[ct][1][{ia[0], ia[1] + 1, 2}] << "|" << std::endl;
		//std::cout << uf_coeff[ct][1][{ia[0], ia[1] + 1, 0}] << ", " << uf_coeff[ct][1][{ia[0], ia[1] + 1, 1}] << ", " << uf_coeff[ct][1][{ia[0], ia[1] + 1, 2}] << "|" << std::endl;

		//std::cout << b[7233] << ", " << (A * p + b)[7233] << std::endl;
		//std::cout << nCell << ", " << nGhost << std::endl;

		//int ia[2] = { ia[0] + 1, ia[1] };
		//auto coeff_pairs = getVelocityUpdateFormula_m<N_DIM>(1, ia, 0, cellType, id, levelSet, ghostId, dc.v_var, jump, sc, ss, dc.blob_material);


		//real err = MLCP_Eigen::MLCP_error(A, flipped_p, b, nCell, nGhost);
		//std::cout << "sym_err " << err << std::endl;

		dc_new.a_var = { getPressure<N>(ss, id, p), ghost_values[0], ghost_values[1], ghost_gaps[0], ghost_gaps[1], id, ghostId[0], ghostId[1] };
		dc_new.levelSet = levelSet0;

		return dc_new;
	}

}
