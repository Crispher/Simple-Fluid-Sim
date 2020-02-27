#pragma once
#include<set>
#include "array.h"
#include "level_set_core.h"

// LESS THAN OR EQUAL TO 0 FOR INSIDE

template<int N>
class LevelSet {
public:
	typedef double real;
	typedef Array<N, real> ArrayND;
private:
	
	ArrayND sdf;

	mutable Array<N, int> *p_id = nullptr;
	
	// optional information
	mutable char flag = 0;

	const int SURFACE = 0b1;
	const int SP = 0b10;
	const int I2SP = 0b100;

	mutable Array<N, int> surface;
	mutable std::vector<std::array<real, N>> sp;
	mutable Array<N, int> i2sp;
	mutable Array<N, int> p_id2;

public:
	LevelSet() {}
	
	LevelSet(const ArrayND&);
	LevelSet(const LevelSet<N>& l) { sdf = l.sdf;}
	void operator = (const LevelSet<N>& l) { sdf = l.sdf;}
	real& operator[] (const int* ia) { return sdf[ia]; }
	real operator[] (const int* ia) const { return sdf[ia]; }

	~LevelSet() { if (p_id != nullptr) delete p_id; };

public:
	Array<N, int> getS() const {
		if (!(flag & SURFACE)) {
			surface = getSDFSignFlip(sdf);
			flag += SURFACE;
		}
		return surface;
	}
	std::vector<std::array<real, N>> getSP() const {
		if (!(flag & SP)) {
			if (!(flag & SURFACE)) {
				surface = getSDFSignFlip(sdf);
				flag += SURFACE;
			}
			sp = initialize_b(sdf, surface);
			flag += SP;
		}
		return sp;
	}
	Array<N, int> getI2SP() const {
		if (!(flag & I2SP)) {
			assert(false);
			std::cout << "use of uninitialized level set" << std::endl;
			exit(0);
		}
		return i2sp;
	}


	//void compute();

	void redistance(real narrowband = 65535.0);

	LevelSet<N> stepMeanCurvatureMotion(real dx, real sigma, real dt, int N) const;

	ArrayND SDF() const {
		return sdf;
	};

	LevelSet<N> setIntersection(const LevelSet<N> &l2) const;
	LevelSet<N> setUnion(const LevelSet<N>& l2) const;
	LevelSet<N> setMinus(const LevelSet<N>& l2) const;
	LevelSet<N> setComplement() const;

	LevelSet<N> operator*(const LevelSet<N> &l2) const { return setIntersection(l2); }
	LevelSet<N> operator+(const LevelSet<N> &l2) const { return setUnion(l2); }
	LevelSet<N> operator-(const LevelSet<N> &l2) const { return setMinus(l2); }
	LevelSet<N> operator!() const { return setComplement(); }

	// l[i,j] is sampled at [x0 + ratio * i, y0 + ratio * j]. resample on [0, 1, ..., dim]
	LevelSet<N> resample(const int *dim, std::vector<real> origin = std::vector<real>(N, 0.0), real ratio = 1) const;

	int identifyComponents() const ;
	// split the level set into connected components
	std::vector<LevelSet<N>> split(int) const;

	bool empty() const { return sdf.min() > 0; };

	ArrayND curvature(real dx) const;

	void save(std::ofstream &os) const;
	void load(std::ifstream &is);

};

template<int N>
inline LevelSet<N>::LevelSet(const ArrayND & s) :
	sdf(s) {
}

template<int N>
inline void LevelSet<N>::redistance(real narrowband) {
	if (!(flag & I2SP)) {
		if (!(flag & SP)) {
			if (!(flag & SURFACE)) {
				surface = getSDFSignFlip(sdf);
				flag += SURFACE;
			}
			sp = initialize_b(sdf, surface);
			flag += SP;
		}
		i2sp = gi2sp(sdf, sp, narrowband);
		flag += I2SP;
	}
}


template<int N>
inline LevelSet<N> LevelSet<N>::setIntersection(const LevelSet<N>& l2) const {
	LevelSet<N> ret = *this;
	for (int i = 0; i < ret.sdf.size(); i++) {
		ret.sdf[i] = std::max<real>(ret.sdf[i], l2.sdf[i]);
	}
	ret.redistance();
	return ret;
}

template<int N>
inline LevelSet<N> LevelSet<N>::setUnion(const LevelSet<N>& l2) const {
	LevelSet<N> ret = *this;
	for (int i = 0; i < ret.sdf.size(); i++) {
		ret.sdf[i] = std::min<real>(ret.sdf[i], l2.sdf[i]);
	}
	ret.redistance();
	return ret;
}

template<int N>
inline LevelSet<N> LevelSet<N>::setMinus(const LevelSet<N>& l2) const {
	LevelSet<N> ret = *this;
	for (int i = 0; i < ret.sdf.size(); i++) {
		ret.sdf[i] = std::max<real>(ret.sdf[i], -l2.sdf[i]);
	}
	ret.redistance();
	return ret;
}

template<int N>
inline LevelSet<N> LevelSet<N>::setComplement() const {
	return LevelSet<N>(sdf * -1.0);
}

template<int N>
inline LevelSet<N> LevelSet<N>::resample(const int * dim, std::vector<real> origin, real ratio) const {
	using namespace range_operator;
	Array<N, real> rs(dim);
	for (auto &[ia, val] : enumerate(rs)) {
		real p[N];
		for (int i = 0; i < N; i++) {
			p[i] = (ia[i] - origin[i]) / ratio;
		}
		val = sdf.interpolate(p);
	}

	return LevelSet<N>(rs * ratio);
}




template<int N>
inline typename LevelSet<N>::ArrayND LevelSet<N>::curvature(real dx) const {
	return getCurvature(sdf) / dx;
}


template<int N>
inline int LevelSet<N>::identifyComponents() const {
	p_id = new Array<N, int>(sdf.dim);
	Array<N, int> &id = *p_id;
	const int UNSEARCHED = -3, INQUEUE = -2, OUTSIDE = -1;
	id = UNSEARCHED;
	int n_component = 0;

	// breadth first search
	using namespace range_operator;
	for (auto[ia, label] : enumerate(id)) {
		if (sdf[ia] > 0) {
			label = OUTSIDE;
			continue;
		}
		if (label == UNSEARCHED) {
			std::queue<int> Q;
			Q.push(ia[N]);
			while (!Q.empty()) {
				int curr = Q.front();
				Q.pop();
				id[curr] = n_component;
				int curr_ia[N];
				id._getIndex(curr, curr_ia);
				auto nbs = id.getNeighborsOf(curr_ia);
				for (int nb : nbs) {
					if (sdf[nb] > 0) {
						id[nb] = OUTSIDE;
						continue;
					}
					if (id[nb] == UNSEARCHED) {
						id[nb] = INQUEUE;
						Q.push(nb);
					}
				}
			}
			n_component++;
		}
		if (label == INQUEUE) {
			assert(false);
		}
	}
	return n_component;
}

template<int N>
inline std::vector<LevelSet<N>> LevelSet<N>::split(int n_component) const {
	std::vector<LevelSet<N>> ret(n_component, LevelSet<N>(sdf));

	using namespace range_operator;
	for (int i = 0; i < n_component; i++) {
		for (auto&[ia, sd] : enumerate(ret[i].sdf)) {
			if ((*p_id)[ia] != i && (*p_id)[ia] >= 0) {
				sd = 10; // clearly outside;
			}
		}
		ret[i].compute();
	}

	return ret;
}


/* implements M. Sussman et.al. 2009: A stable and efficient method
for treating surface tension in incompressible two-phase flow */
template<int N>
inline LevelSet<N> LevelSet<N>::stepMeanCurvatureMotion(real dx, real sigma, real dt, int N) const {
	using namespace range_operator;

	int sevens[N], threes[N];
	for (int j = 0; j < N; j++) {
		sevens[j] = 7;
		threes[j] = 3;
	}
	Array<N, real> p_weight(sevens);
	for (auto&[ia, w] : enumerate(p_weight)) {
		real dsq = 0;
		for (int j = 0; j < N; j++) {
			dsq += (ia[j] - 3) * (ia[j] - 3);
		}
		w = 1.0 / pow(dsq, 5);
	}
	p_weight[threes] = 0.0;

	Array<N, real> SUM(sdf.dim), W(sdf.dim), S = sdf;

	auto surface = getSDFSignFlip(S);

	for (int i = 0; i < N; i++) {
		
		// auto surface = getSDFSignFlip(S); // s = S * dx
		auto Kappa = getCurvature(S, surface); // kappa = Kappa / dx;
		
		// real kappa_avg = averageCurvature(kappa, s);
		real Kappa_avg = (Kappa * surface).sum() / surface.sum();

		// extrapolate kappa
		W = 0; SUM = 0;
		auto offset = dim2offset<N>(S.dim);

		for (int j = 0; j < Kappa.size(); j++) {
			if (!surface[j])
				continue;
			auto ia = i2ia<N>(j, offset);
			real k = Kappa[j];
			for (auto local_ia : Range<N>(sevens)) {
				int ref_ia[N];
				bool valid = true;
				for (int j = 0; j < N; j++) {
					ref_ia[j] = ia[j] + local_ia[j] - 3;
					if (ref_ia[j] < 0 || ref_ia[j] >= sdf.dim[j]) {
						valid = false;
					}
				}
				if (!valid)
					continue;
				W[ref_ia] += p_weight[local_ia];
				SUM[ref_ia] += p_weight[local_ia] * k;
			}
		}


		for (int j = 0; j < Kappa.size(); j++) {

			if (W[j] > 0 && !surface[j]) {
				Kappa[j] = SUM[j] / W[j];
			}
		}

		// evolve d
		Array<N, real> D = (Kappa - Kappa_avg) * (sigma * dt / dx / dx);
		S += D;
		auto abs = [](real r) -> real {return r > 0 ? r : -r; };
		real Dmax = D.max(abs);

		// reinitialize d, Sect. 6

		// step 1.
		// LevelSet<N_DIM> sl = S;
		surface = getSDFSignFlip(S);
		if (i < N - 1)
			fastMarching<N, real>(S, std::max<real>(Dmax + 1, 2), surface);
		else
			fastMarching<N, real>(S, 10, surface);
		// S = sl.sdf;
	}
	// getchar();

	LevelSet<N> ret(S);
	return ret;
}

template<int N>
inline void LevelSet<N>::save(std::ofstream & os) const {
	os << "levelset-begin" << std::endl;
	sdf.save(os);
	os << "levelset-end" << std::endl;
}

template<int N>
inline void LevelSet<N>::load(std::ifstream & is) {
	std::string s;
	is >> s;
	if (s != "levelset-begin") {
		std::cout << "failed to load level set: " << s << std::endl;
	}
	sdf.load(is);
	is >> s;
	if (s != "levelset-end") {
		std::cout << "failed to load level set: " << s << std::endl;
	}
}



#pragma region LEVEL_SET_UTILITIES
template<int N, typename real>
Array<N, real> getWall(const int *dim) {
	Array<N, LevelSet<N>::real> sdf(dim);
	sdf = 0.5;
	using namespace range_operator;
	for (int i = 0; i < N; i++) {
		sdf[R<N>(i, { 0 })] = -0.5;
		sdf[R<N>(i, { -1 })] = -0.5;
	}
	fastMarching(sdf, 10.0);
	return sdf;
}

template<int N, typename real>
Array<N, real> getSphere(const int *dim, const real *center, typename LevelSet<N>::real radius) {
	Array<N, LevelSet<N>::real> sdf([&](const int *ia) {
		LevelSet<N>::real d = 0;
		for (int i = 0; i < N; i++) {
			d += (ia[i] - center[i]) * (ia[i] - center[i]);
		}
		return sqrt(d) - radius;
	}, dim);

	return sdf;
}

template<int N, typename real>
Array<N, real> getCuboid(const int *dim, const real *left, const real *right) {
	Array<N, typename LevelSet<N>::real> l(dim);
	l = -1;
	LevelSet<N> ret = l;
	for (int i = 0; i < N; i++) {
		LevelSet<N> tl = Array<N, typename LevelSet<N>::real>([&](const int *ia) { return left[i] - ia[i]; }, dim);
		ret = ret * tl;
		ret = ret * LevelSet<N>(Array<N, typename LevelSet<N>::real>([&](const int *ia) {return ia[i] - right[i]; }, dim));
	}
	return ret.SDF();
}


template<int N, typename R>
Array<N, int> getSignFlip(const Array<N, R> &s) {
	using namespace range_operator;

	Array<N, int> ret(s.dim, 0);
	
	for (auto ia : Range<N>(s.dim)) {
		for (int i = 0; i < N; i++) {
			if (ia[i] + 1 < s.dim[i]) {
				R sl = s[ia];
				R sr = s[I<N>(ia, i, 1)];
				if (!same_sign(sl, sr)) {
					ret[ia] = 1;
					ret[I<N>(ia, i, 1)] = 1;
				}
			}
		}
	}
	return ret;
}

template<int N, typename R>
Array<N, int> getSDFSignFlip(const Array<N, R> &s, R clamp = 2) {
	using namespace range_operator;

	Array<N, int> ret(s.dim);
	ret = 0;
	for (auto ia : Range<N>(s.dim)) {
		if (abs(s[ia]) > clamp)
			continue;
		for (int i = 0; i < N; i++) {
			if (ia[i] + 1 < s.dim[i]) {
				R sl = s[ia];
				R sr = s[I<N>(ia, i, 1)];
				if (!same_sign(sl, sr)) {
					ret[ia] = 1;
					ret[I<N>(ia, i, 1)] = 1;
				}
			}
		}
	}
	return ret;
}


template<typename real>
inline Array<2, real> getCurvature(const Array<2, real> &sdf) {
	return getCurvature(sdf, getSignFlip(sdf));
}

template<typename real>
inline Array<3, real> getCurvature(const Array<3, real> &sdf) {
	return getCurvature(sdf, getSDFSignFlip(sdf));
}

template<int N, typename real>
inline void fastMarching(Array<N, real> &sdf, real narrow_band) {
	fastMarching<N, real>(sdf, narrow_band, getSDFSignFlip(sdf));
}




template<int N, typename real>
void extrapolate(Array<N, int> &valid, Array<N, real> &u, int n_sweeps) {
	using namespace range_operator;
	std::vector<int> Q;
	auto offset = dim2offset<N>(u.dim);

	for (auto ia : Range<N>(u.dim)) {
		int id = ia2i<N>(ia, offset);
		if (valid[ia] == 0) {
			for (auto nb : u.getNeighborsOf(ia)) {
				if (valid[nb] == 1) {
					Q.push_back(id);
					valid[ia] = 2;
					break;
				}
			}
		}
	}

	
	for (int L = 0; L < n_sweeps; L++) {
		
		if (Q.empty()) {
			break;
		}
		std::vector<int> Q2;
		for (auto id : Q) {
			auto ia = i2ia<N>(id, offset);
			real S = 0.0, W = 0.0;
			for (auto nb : u.getNeighborsOf(ia)) {
				if (valid[nb] == 1) {
					S += u[nb];
					W += 1.0;
				} else if (valid[nb] == 0) {
					valid[nb] = 2;
					Q2.push_back(nb);
				}
			}
			if (W == 0) std::cout << "EEEEEEEEEEEEEEEEEEEEEESRRRRRRRRRRRRRRRRRRRRRRR" << std::endl;
			u[id] = S / W;
		}
		for (auto id : Q) {
			valid[id] = 1;
		}
		Q = Q2;
	}

}

template<int N, typename real>
Array<N, real> advect(const Array<N, real> &sdf, const std::vector<std::array<real, N>> &sp, const std::vector<std::array<real, N>> &v_sp, const Array<N, int> &i2sp, real dt, real dx) {
	using namespace range_operator;
	auto offset = dim2offset<N>(sdf.dim);

	auto ret = sdf;

	for (int i = 0; i < sdf.size(); i++) {
		if (i2sp[i] < 0)
			continue;
		auto ia = i2ia<N>(i, offset);
		real tb[N];
		for (int j = 0; j < N; j++) {
			tb[j] = ia[j] - v_sp[i2sp[i]][j] * dt / dx;
		}
		ret[i] = sdf.interpolate_cubic<decltype(tb), false>(tb);
	}

	return ret;
}

// incorrect yet
template<int N, typename real>
std::vector<std::array<real, N>> initialize(Array<N, real> &sdf, const Array<N, int> &surface) {
	using namespace range_operator;
	const real d_th = 1e-5;
	auto offset = dim2offset<N>(sdf.dim);

	Array<N, real> normal[N];
	for (int i = 0; i < N; i++) {
		normal[i] = Array<N, real>(sdf.dim);
		for (auto &[ia, val] : enumerate(normal[i])) {
			if (ia[i] == 0) {
				val = sdf[I<N>(ia, i, 1)] - sdf[ia];
			} else if (ia[i] + 1 == sdf.dim[i]) {
				val = sdf[ia] - sdf[I<N>(ia, i, -1)];
			} else {
				val = (sdf[I<N>(ia, i, 1)] - sdf[I<N>(ia, i, -1)]) / 2.0;
			}
		}
	}


	std::vector<std::array<real, N>> ret;
	std::vector<real> sdf_update;
	for (int id = 0; id < surface.size(); id++) {
		if (surface[id]) {
			auto ia = i2ia<N>(id, offset);
			// real pos[N];
			std::array<real, N> pos;
			copy_indexable<N>(ia, pos);
			
			real n[N], n_sq = 0;
			int c[N];
			real d = sdf[id];
			if (abs(d) <= d_th) {
				sdf_update.push_back(d);
				ret.push_back(pos);
				continue;
			}
			int cnt = 0;
			while (abs(d) > d_th) {
				for (int i = 0; i < N; i++) {
					c[i] = std::clamp<int>((int)(pos[i]), 0, sdf.dim[i] - 2);
				}
				for (int i = 0; i < N; i++) {
					n[i] = normal[i].interpolate(pos);
					n_sq += n[i] * n[i];
				}
				real g_norm = sqrt(n_sq);
				for (int i = 0; i < N; i++) {
					pos[i] -= d * n[i] / g_norm;
				}
				d = sdf.interpolate(pos);
				cnt++;
				if (cnt > 10000) {
					std::cout << "LL" << pos[0] << ", " << pos[1] << ": " << d << std::endl;
				}
			}
			real sd = 0;
			for (int i = 0; i < N; i++) {
				sd += (pos[i] - (real)ia[i]) * (pos[i] - (real)ia[i]);
			}
			sd = sdf[id] >= 0 ? sqrt(sd) : -sqrt(sd);
			sdf_update.push_back(sd);
			ret.push_back(pos);
			
		}
	}
	int cursor = 0;
	for (int id = 0; id < sdf.size(); id++) {
		if (surface[id]) {
			sdf[id] = sdf_update[cursor++];
		}
	}
	return ret;
}

template<int N, typename real>
Array<N, int> march(Array<N, real> &sdf, const Array<N, int> &surface, std::vector<std::array<real, N>> sp) {
	using namespace range_operator;
	
	Array<N, real> usdf(sdf.dim, 2.0);
	Array<N, int> i2sp(sdf.dim, -1);
	int cursor = 0;
	for (int i = 0; i < usdf.size(); i++) {
		if (surface[i]) {
			usdf[i] = abs(sdf[i]);
			i2sp[i] = cursor++;
		}
	}


	// construct queue
	auto comparator = [&usdf](int a, int b) -> bool {
		if ((usdf[a]) < (usdf[b])) {
			return true;
		}
		if ((usdf[a]) > (usdf[b])) {
			return false;
		}
		return a < b;
	};
	std::set<int, decltype(comparator)> Q(comparator);
	int const UNKNOWN = -2, const INQUEUE = -1, const KNOWN = 0;
	Array<N, int> status(sdf.dim, UNKNOWN);



	for (int i = 0; i < usdf.size(); i++) {
		if (surface[i]) {
			status[i] = INQUEUE;
			Q.insert(i);
		}
	}

	// iterate through all grid points
	real max_d = 0;
	// std::vector<int> neighbors;
	while (!Q.empty() && max_d < 10) {
		int p = *(Q.begin());
		status[p] = KNOWN;
		real sign = (sdf[p] >= 0) ? 1.0 : -1.0;
		sdf[p] = sign * abs(usdf[p]);
		Q.erase(p);
		max_d = usdf[p];
		int pia[N];
		usdf._getIndex(p, pia);
		auto neighbors = usdf.getNeighborsOf_O2(pia);
		for (int k = 0; true; k++) {
			int nb = neighbors[k];
			if (nb < 0) break;
			real dsq = 0.0;
			int nbia[N];
			status._getIndex(nb, nbia);
			for (int i = 0; i < N; i++) {
				dsq += (nbia[i] - sp[i2sp[p]][i]) * (nbia[i] - sp[i2sp[p]][i]);
			}
			switch (status[nb]) {
			case UNKNOWN:
				// push to queue
				usdf[nb] = sqrt(dsq);
				status[nb] = INQUEUE;
				Q.insert(nb);
				i2sp[nb] = i2sp[p];
				break;
			case INQUEUE:
				// update if necessary
				if (dsq < usdf[nb] * usdf[nb]) {
					Q.erase(nb);
					usdf[nb] = sqrt(dsq);
					Q.insert(nb);
					i2sp[nb] = i2sp[p];
				}
				break;
			case KNOWN:
				break;
			}
		}
	}
	for (auto[ia, s] : enumerate(status)) {
		if (s != KNOWN) {
			sdf[ia] = sdf[ia] > 0 ? 10.0 : -10.0;
		}
	}
	return i2sp;
}
#pragma endregion