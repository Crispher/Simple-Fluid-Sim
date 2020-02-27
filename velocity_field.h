#pragma once
#include "array.h"
#include<iomanip>
#include "opencv2/opencv.hpp"
#include "particle.h"

template <int N>
class LevelSet;

/* The field is defined on [-0.5, dim[i]-0.5], while individual component is defined on [-0.5, dim[i]-0.5]
v values are w.r.t. real world dimension (not integer grid).
*/

template <int N>
class VelocityField {
public:
	typedef double real;
	typedef Array<N, real> ArrayND;
	int dim[N];
	// private:

	ArrayND v[N];

	real dx;

public:
	template<typename Indexable>
	VelocityField(Indexable dims, real dx);
	VelocityField() {}
	VelocityField(std::initializer_list<int> dims, real dx);
	VelocityField(std::ifstream &is) { load(is); }

	VelocityField(const VelocityField<N>&);
	void operator=(const VelocityField<N>&);

public:
	// advect phi, where phi[i,j] is sampled at (x0 + i * ratio, y0 + j * ratio), (x0, y0) being 'origin'.
	ArrayND advect(const ArrayND&, real dt, std::vector<real> origin = std::vector<real>(N, 0.0), real ratio = 1.0, const Array<N, real> *v_at_sample = nullptr) const;

	// advect phi, where phi[i,j] is sampled at (x0 + i * ratio, y0 + j * ratio), (x0, y0) being 'origin'.
	template<bool clamp = true>
	ArrayND advect_cubic(const ArrayND&, real dt, std::vector<real> origin = std::vector<real>(N, 0.0), real ratio = 1.0) const;

	Array<N, real> advectLevelSet(const Array<N, real>& phi, real dt, real narrow_band, const Array<N, int> &i2sp, const std::vector<std::array<real, N>> &v_sp);
	
	LevelSet<N> advectLevelSet(const LevelSet<N> &phi, real dt, real narrow_band);

	Array<N, std::array<real, N>> traceback_cubic(real dt) const;

	VelocityField<N> advectVelocity(const VelocityField<N>&, real dt) const;
	VelocityField<N> advectVelocity_cubic(const VelocityField<N>&, real dt) const;

	void getVelocity(const std::vector<Particle<N>> & particles);
	void getParticleVelocity(const std::vector<Particle<N>> & particles);
	void advectParticle(std::vector<Particle<N>> & particles, real dt) const;
	void transferVelocity(std::vector<Particle<N>> & particles) const;
	void incrementParticleVelocity(std::vector<Particle<N>> &particles) const;
	
	VelocityField<N> extrapolate_v(const Array<N, real> &sdf) const;

	template<typename real_vec>
	std::array<real, N> interpolate(const real_vec &pos) const;
	template<typename real_vec>
	std::array<real, N> interpolate_cubic(const real_vec &pos) const;

	ArrayND divergence() const;

	std::array<Array<N, real>, N> getCellCenterAvg() const;

	void save(std::ofstream &os) const;
	void load(std::ifstream &os);

	real getCFLTimeStep(real CFL_const = 1.0) const;
	ArrayND& operator[](int i) {
		return v[i];
	}
	const ArrayND& operator[](int i) const {
		return v[i];
	};

public:
	ArrayND applyViscosityComponent(int axis, real viscosity, real dt, const Array<N, int>&);
	VelocityField<N> applyViscosity(real viscosity, real dt, const Array<N, int>& ct);

public:
	void operator+=(Array<1, real>);


	void setZeroBoundary();

	template<typename real_vec>
	std::vector<real> _getInterpolatePos(int C, const real_vec &pos) const;
};



template<int N>
template<typename Indexable>
inline VelocityField<N>::VelocityField(Indexable dims, real _dx) {
	for (int i = 0; i < N; i++) {
		//dim[i] = dims[i];
		int d = dims[i];
		dim[i] = d;
	}
	for (int i = 0; i < N; i++) {
		dim[i]++;
		v[i] = ArrayND(dim);
		dim[i]--;
	}
	dx = _dx;
}

template<int N>
inline VelocityField<N>::VelocityField(std::initializer_list<int> dims, real _dx) :
	VelocityField(std::vector<int>(dims), _dx) 		{
}

template<int N>
inline VelocityField<N>::VelocityField(const VelocityField<N>& _v) {
	for (int i = 0; i < N; i++) {
		dim[i] = _v.dim[i];
		v[i] = _v.v[i];
	}
	dx = _v.dx;
}

template<int N>
inline void VelocityField<N>::operator=(const VelocityField<N>& _v) {
	//todo: assert dimension agree
	for (int i = 0; i < N; i++) {
		dim[i] = _v.dim[i];
		v[i] = _v.v[i];
	}
	dx = _v.dx;
}

template<int N>
inline typename VelocityField<N>::ArrayND VelocityField<N>::advect(const ArrayND &phi, real dt, std::vector<real> origin, real ratio, const Array<N, real> *v_at_sample) const {
	using namespace range_operator;
	ArrayND advected(phi.dim);

	// optimization
	bool optimize = true;
	for (int i = 0; i < N; i++) {
		optimize = optimize && (origin[i] == 0);
	}
	optimize = optimize && (ratio == 1);
	Array<N, real> vas[N];
	for (int i = 0; i < N; i++) {
		if (optimize) {
			vas[i] = (v[i][R<N>(i, { _, -1 })] + v[i][R<N>(i, { 1, _ })]) / 2;
			v_at_sample = vas;
		} else {
			vas[i].disableWarning();
		}
	}

	std::vector<real> interpolatePos(N);
	std::array<real, N> velocity;
	for (auto &[ia, phi_advected] : enumerate(advected)) {
		// compute new value for phi[ia]	
		for (int i = 0; i < N; i++) {
			interpolatePos[i] = origin[i] + ia[i] * ratio;
		}
		if (v_at_sample != nullptr) {
			// for optimization purpose
			for (int i = 0; i < N; i++) {
				velocity[i] = v_at_sample[i][ia];
			}
		} else {
			velocity = interpolate(interpolatePos);
		}
		for (int i = 0; i < N; i++) {
			interpolatePos[i] -= velocity[i] * dt / dx;
			interpolatePos[i] -= origin[i];
			interpolatePos[i] /= ratio;
		}
		phi_advected = phi.interpolate(interpolatePos);
	}
	return advected;
}

template<int N>
template<bool clamp>
inline typename VelocityField<N>::ArrayND VelocityField<N>::advect_cubic(const ArrayND &phi, real dt, std::vector<real> origin, real ratio) const {
	using namespace range_operator;
	ArrayND advected(phi.dim);

	std::array<real, N> interpolatePos;
	std::array<real, N> velocity;
	for (auto &[ia, phi_advected] : enumerate(advected)) {
		if (abs(phi_advected) > 6)
			continue;
		// compute new value for phi[ia]	
		for (int i = 0; i < N; i++) {
			interpolatePos[i] = origin[i] + ia[i] * ratio;
		}
		velocity = interpolate_cubic(interpolatePos);
		for (int i = 0; i < N; i++) {
			interpolatePos[i] -= velocity[i] * dt / dx;
			interpolatePos[i] -= origin[i];
			interpolatePos[i] /= ratio;
		}
		phi_advected = phi.interpolate_cubic<decltype(interpolatePos), clamp>(interpolatePos);
	}
	return advected;
}

template<int N>
inline typename VelocityField<N>::ArrayND VelocityField<N>::advectLevelSet(const Array<N, real>& phi, real dt, 
		real narrow_band, const Array<N, int>& i2sp, const std::vector<std::array<real, N>>& v_sp) {
	using namespace range_operator;
	auto ret = phi;
	std::array<real, N> interp_pos;
	auto offset = dim2offset<N>(phi.dim);

	for (int i = 0; i < phi.size(); i++) {
		if (abs(phi[i]) >= narrow_band) {
			continue;
		}
		auto ia = i2ia<N>(i, offset);
		auto &v = v_sp[i2sp[i]];
		for (int j = 0; j < N; j++) {
			interp_pos[j] = ia[j] - v[j] * dt / dx;
		}
		ret[i] = phi.interpolate(interp_pos);
	}
	return ret;
}

template<int N>
inline LevelSet<N> VelocityField<N>::advectLevelSet(const LevelSet<N>& phi, real dt, real narrow_band) {
	using namespace range_operator;

	auto sp = phi.getSP();
	std::vector<std::array<real, N>> v_sp(sp.size());
	for (int i = 0; i < sp.size(); i++) {
		v_sp[i] = interpolate(sp[i]);
	}

	auto sdf = advectLevelSet(phi.SDF(), dt, narrow_band, phi.getI2SP(), v_sp);

	LevelSet<N> ret(sdf);
	ret.redistance();
	return ret;
}

template<int N>
inline Array<N, std::array<typename VelocityField<N>::real, N>> VelocityField<N>::traceback_cubic(real dt) const {
	using namespace range_operator;
	Array<N, std::array<real, N>> tb(dim);

	// optimization
	Array<N, real> vas[N];
	for (int i = 0; i < N; i++) {
		vas[i] = (v[i][R<N>(i, { _, -1 })] + v[i][R<N>(i, { 1, _ })]) / 2;
	}
	// end optimization

	std::array<real, N> velocity;
	for (auto &[ia, tbp] : enumerate(tb)) {

		for (int i = 0; i < N; i++) {
			velocity[i] = vas[i][ia];
		}

		for (int i = 0; i < N; i++) {
			tbp[i] -= velocity[i] * dt / dx;
		}
	}
	return tb;
}

template<int N>
inline VelocityField<N> VelocityField<N>::advectVelocity(const VelocityField<N> &_v, real dt) const {
	VelocityField<N> ret(dim, dx);
	using namespace range_operator;
	for (int i = 0; i < N; i++) {
		std::vector<real> origin(N, 0.0);
		origin[i] = -0.5;
		Array<N, real> velocity[N];
		for (int j = 0; j < N; j++) {
			velocity[j] = Array<N, real>(_v[i].dim);
			if (j == i) {
				velocity[j] = v[j];
			} else {
				Array<N, real> vj = ((v[j][R<N>(i, { _, -1 }, j, { _, -1 })] + v[j][R<N>(i, { 1, _ }, j, { _, -1 })]
					+ v[j][R<N>(i, { _, -1 }, j, { 1, _ })] + v[j][R<N>(i, { 1, _ }, j, { 1, _ })]) * 0.25);
				velocity[j] = vj.pad(i, 1, 1, array_util::COPY_OUTERMOST);
			}
		}
		ret[i] = advect(_v[i], dt, origin, 1.0, velocity);
	}
	return ret;
}

template<int N>
inline VelocityField<N> VelocityField<N>::advectVelocity_cubic(const VelocityField<N>& _v, real dt) const {
	VelocityField<N> ret(dim, dx);
	using namespace range_operator;
	for (int i = 0; i < N; i++) {
		std::vector<real> origin(N, 0.0);
		origin[i] = -0.5;
		ret[i] = advect_cubic(_v[i], dt, origin, 1.0);
	}
	return ret;
}

template<int N>
inline VelocityField<N> VelocityField<N>::extrapolate_v(const Array<N, real>& sdf) const {
	VelocityField<N> ret = *this;

	Array<N, int> valid[N];
	auto le0 = [](real r) -> int { return r <= 0; };
	Array<N, int> cValid = sdf.map<int>(le0);

	auto binarize = [](int c) -> int { return c > 0; };
	for (int i = 0; i < N; i++) {
		valid[i] = (cValid.pad(i, 0, 1, true) + cValid.pad(i, 1, 0, true)).map<int>(binarize);
	}

	for (int i = 0; i < N; i++) {
		extrapolate(valid[i], ret.v[i], 30);
	}

	return ret;
}

template<int N>
template<typename real_vec>
inline std::array<typename VelocityField<N>::real, N> VelocityField<N>::interpolate(const real_vec &pos) const {
	std::array<real, N> ret;
	for (int i = 0; i < N; i++) {
		ret[i] = v[i].interpolate(_getInterpolatePos(i, pos));
	}
	return ret;
}

template<int N>
template<typename real_vec>
inline std::array<typename VelocityField<N>::real, N> VelocityField<N>::interpolate_cubic(const real_vec & pos) const {
	std::array<real, N> ret;
	for (int i = 0; i < N; i++) {
		ret[i] = v[i].interpolate_cubic(_getInterpolatePos(i, pos));
	}
	return ret;
}

template<int N>
inline typename VelocityField<N>::ArrayND VelocityField<N>::divergence() const {
	using underscore_as_semicolon::_;
	ArrayND div(dim);
	div = 0.0;
	for (int i = 0; i < N; i++) {
		auto l = get_full_range<N>(), r = l;
		l[i] = { 0, -1 }; r[i] = { 1, _ };
		div += (v[i][r] - v[i][l]) / dx;
	}
	return div;
}

template<int N>
inline std::array<Array<N, typename VelocityField<N>::real>, N> VelocityField<N>::getCellCenterAvg() const {
	using namespace range_operator;
	using underscore_as_semicolon::_;
	std::array<Array<N, real>, N> ret;

	for (int i = 0; i < N; i++) {
		ret[i] = (v[i][R<N>(i, { _, -1 })] + v[i][R<N>(i, { 1, _ })]) / 2.0;
	}
	return ret;
}

template<int N>
inline void VelocityField<N>::save(std::ofstream & os) const {
	os << "velocityfield-begin" << std::endl;
	for (int i = 0; i < N; i++) {
		os << dim[i] << std::endl;
	}
	os << dx << std::endl;
	for (int i = 0; i < N; i++) {
		v[i].save(os);
	}
	os << "velocityfield-end" << std::endl;
}

template<int N>
inline void VelocityField<N>::load(std::ifstream & is) {
	std::string s;
	is >> s;
	if (s == "velocityfield-begin") {
		for (int i = 0; i < N; i++) {
			is >> dim[i];
		}
		is >> dx;
		for (int i = 0; i < N; i++) {
			v[i].load(is);
		}
	}
	is >> s;
	if (s != "velocityfield-end") {
		std::cout << "failed to load velocity field" << std::endl;
	}
}

template<int N>
inline typename VelocityField<N>::real VelocityField<N>::getCFLTimeStep(real CFL_const) const {
	real vm = 1e-5;
	for (int i = 0; i < N; i++) {
		vm = std::max<real>(vm, v[i].max());
		vm = std::max<real>(vm, -v[i].min());
	}
	return CFL_const / vm;
}

template<int N>
inline VelocityField<N> VelocityField<N>::applyViscosity(real viscosity, real dt, const Array<N, int>& ct) {
	VelocityField<N> ret(dim, dx);
	for (int i = 0; i < N; i++) {
		ret.v[i] = applyViscosityComponent(i, viscosity, dt, ct);
	}
	return ret;
}

template<int N>
inline void VelocityField<N>::setZeroBoundary() {
	using namespace range_operator;
	for (int i = 0; i < N; i++) {
		v[i][R<N>(i, { 0 })] = 0;
		v[i][R<N>(i, { -1 })] = 0;
	}
}

template<int N>
template<typename real_vec>
inline std::vector<typename VelocityField<N>::real> VelocityField<N>::_getInterpolatePos(int C, const real_vec& pos) const {
	std::vector<real> ret(N);
	copy_indexable<N>(pos, ret);
	ret[C] += 0.5;
	return ret;
}

template<int N>
inline void VelocityField<N>::getVelocity(const std::vector<Particle<N>>& particles) {
	using range_operator::R;

	int dim1[N];
	for (int i = 0; i < N; i++) {
		dim1[i] = dim[i] + 1;
	}
	Array<N, std::vector<int>> nearbyParticles(dim1);
	ArrayND sdf_sign(dim);

	for (auto p = particles.begin(); p != particles.end(); p++) {
		// add to nearby particles
		// assume particles have legal positions (within the domain)
		auto nearbyCells = R<N>();
		for (int i = 0; i < N; i++) {
			int l = p->pos[i] - 0.5 - 1e-4;
			int r = p->pos[i] + 1.5 + 1e-4;
			l = std::max<int>(l, 0);
			r = std::min<int>(r + 1, dim[i]);
			nearbyCells[i] = { l, r };
		}
		nearbyParticles[nearbyCells].forEach([&](int *ia, std::vector<int> &v) {
			v.push_back(p - particles.begin());
		});

	}

	auto positionalW = [](real* gp, const std::array<real, N> &pp, real radius) -> real {
		real ret = 1.0;
		for (int i = 0; i < N; i++) {
			ret *= radius - std::abs<real>(gp[i] - pp[i]);
			if (ret <= 0) return 0;
		}

		return ret;
	};

	sdf_sign.forEach([&](int *ia, real &val) {
		std::vector<int> ps = nearbyParticles[ia];
		val = 1.0;
		for (int j : ps) {
			real gridpoint[N];
			copy_indexable<N>(ia, gridpoint);
			if (positionalW(gridpoint, particles[j].pos, 0.5 + 1e-3) > 0) {
				val = -1.0;
			}
		}
	});

	for (int i = 0; i < N; i++) {
		v[i].forEach([&](int *ia, real &val) {
			std::vector<int> ps = nearbyParticles[ia];
			if (!ps.empty()) {
				real gridpoint[N];
				copy_indexable<N>(ia, gridpoint);
				gridpoint[i] -= 0.5;
				std::vector<real> w(ps.size(), 1);
				real w_sum = 0;
				val = 0;
				for (int j = 0; j < ps.size(); j++) {
					w[j] = positionalW(gridpoint, particles[ps[j]].pos, 1.0 + 1e-3);
					w_sum += w[j];
					val += particles[ps[j]].vel[i] * w[j];
				}
				if (w_sum < 1e-5) {
					std::cout << "0 weight " << w_sum << ", " << sdf_sign[ia] << std::endl;
					std::cout << w[0] << ", " << w[1] << std::endl;
				} else {
					val /= w_sum;
				}
			}
		});
	}

	print();
	sdf_sign.print();
	// extrapolate

	// auto _v = extrapolateVelocity(sdf_sign);
	print();
}

template<int N>
inline void VelocityField<N>::getParticleVelocity(const std::vector<Particle<N>>& particles) {
	using namespace range_operator;
	int dim1[N];
	for (int i = 0; i < N; i++) {
		dim1[i] = dim[i] + 1;
	}
	Array<N, std::vector<int>> p1 = get_nearby_particles(dim1, 1.5 + 1e-3, particles);
	Array<N, std::vector<int>> p2 = get_nearby_particles(dim1, 6.5, particles);

	auto positionalW = [](real* gp, const std::array<real, N> &pp, real radius) -> real {
		real ret = 1.0;
		for (int i = 0; i < N; i++) {
			ret *= radius - std::abs<real>(gp[i] - pp[i]);
			if (ret <= 0) return 0;
		}
		return ret;
	};

	auto distance = [](real *gp, const std::array<real, N> &pp) -> real {
		real ret = 0;
		for (int i = 0; i < N; i++) {
			ret += (gp[i] - pp[i]) * (gp[i] - pp[i]);
		}
		return ret;
	};

	for (int i = 0; i < N; i++) {
		// std::cout << i << std::endl;
		for (auto &[ia, val] : enumerate(v[i])) {
			real gridpoint[N];
			copy_indexable<N>(ia, gridpoint);
			gridpoint[i] -= 0.5;

			std::vector<int> r1 = p1[ia];
			real v_sum = 0;
			real w_sum = 0;
			for (int id : r1) {
				// std::cout << id << std::endl;
				real wid = positionalW(gridpoint, particles[id].pos, 1.0);
				v_sum += wid * particles[id].vel[i];
				w_sum += wid;
			}
			if (w_sum > 1e-5) {
				// if there are valid particles
				val = v_sum / w_sum;
			} else {
				// if no particles within r=1, use vel of the nearest particle
				std::vector<int> r2 = p2[ia];
				real dmax = 10000;
				for (int id : r2) {
					// std::cout << id << "-" << std::endl;
					real d = distance(gridpoint, particles[id].pos);
					if (d < dmax) {
						val = particles[id].vel[i];
						dmax = d;
					}
				}
			}
		}
	}

	// print();
}

template<int N>
inline void VelocityField<N>::advectParticle(std::vector<Particle<N>> & particles, real dt) const {
	for (auto &p : particles) {
		// advect p, RK2
		//std::vector<real> vp = interpolate(p.pos);
		std::array<real, N> vp = interpolate(p.pos);
		//std::vector<real> pp = p.pos;
		std::array<real, N> pp = p.pos;
		for (int i = 0; i < N; i++) {
			pp[i] += vp[i] * dt / dx / 2.0;
			pp[i] = std::clamp<real>(pp[i], -0.5 + 1e-3, dim[i] - 0.5 + 1e-3);
		}
		std::array<real, N> vp2 = interpolate(pp);
		for (int i = 0; i < N; i++) {
			p.pos[i] += vp2[i] * dt / dx;
			p.pos[i] = std::clamp<real>(p.pos[i], -0.5 + 1e-3, dim[i] - 0.5 + 1e-3);
		}
	}

}

template<int N>
inline void VelocityField<N>::transferVelocity(std::vector<Particle<N>>& particles) const {
	for (auto &p : particles) {
		// interpolate to p.pos
		// std::array<real, N_DIM> pp;
		// p.vel = interpolate(p.pos);
		p.vel = interpolate(p.pos);
	}
}

template<int N>
inline void VelocityField<N>::incrementParticleVelocity(std::vector<Particle<N>>& particles) const {
	for (auto &p : particles) {
		// auto DV = interpolate(p.pos);//todo
		auto DV = interpolate(p.pos);
		for (int i = 0; i < N; i++)
			p.vel[i] += DV[i];
	}
}
