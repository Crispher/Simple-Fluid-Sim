#pragma once
#include "array.h"

inline bool same_sign(double l, double r) {
	return (l > 0) && (r > 0) || (l <= 0) && (r <= 0);
}

inline double length_fraction(double l, double r) {
	if (l == 0 && r == 0) {
		return 0.5;
	}
	if (l <= 0 && r <= 0) {
		return 1.0;
	}
	if (l >= 0 && r >= 0) {
		return 0.0;
	}
	if (l < 0 && r > 0) {
		return l / (l - r);
	}
	if (l > 0 && r < 0) {
		return r / (r - l);
	}
}

template<typename real>
inline real hf_curv2d(real l, real m, real r) {
	return (l - 2 * m + r) / sqrt(1 + ((r - l) / 2)*((r - l) / 2));
}

template<typename real>
inline real hf_curv3d(real H[3][3], real gamma = 0.0) {
	const int l = 0, m = 1, r = 2;
	real Hx = (gamma * (H[r][r] - H[l][r]) + H[l][m] - H[l][m] + gamma * (H[r][l] - H[l][l])) / 2.0 / (1 + 2 * gamma);
	real Hxx = (gamma * (H[r][r] - 2 * H[m][r] + H[l][r]) + H[r][m] - 2 * H[m][m] + H[l][m] +
		gamma * (H[r][l] - 2 * H[m][l] + H[l][l])) / (1 + 2 * gamma);
	real Hxy = (H[r][r] - H[l][r] - H[r][l] + H[l][l]) / 4.0;
	real Hy = (gamma * (H[r][r] - H[r][l]) + H[m][l] - H[m][l] + gamma * (H[l][r] - H[l][l])) / 2.0 / (1 + 2 * gamma);
	real Hyy = (gamma * (H[r][r] - 2 * H[r][m] + H[r][l]) + H[m][r] - 2 * H[m][m] + H[m][l] +
		gamma * (H[l][r] - 2 * H[l][m] + H[l][l])) / (1 + 2 * gamma);
	real kappa = (Hxx + Hyy + Hxx * Hy * Hy + Hyy * Hx * Hx - 2 * Hxy * Hx * Hy) / pow(1 + Hx * Hx + Hy * Hy, 1.5);
	return kappa;
}

// return the curvature of a 2D level set, for cells near the surface.
template<typename real>
inline Array<2, real> getCurvature(const Array<2, real> &sdf, const Array<2, int> &surface) {
	using namespace range_operator;
	Array<2, real> ret(sdf.dim, 0.0);

	int S = 2;
	for (auto ia : Subrange<2>(sdf.dim, { { S, -S },{ S, -S } })) {
		if (!surface[ia])
			continue;

		real DX = sdf(ia[0] + 1, ia[1]) - sdf(ia[0] - 1, ia[1]);
		real DY = sdf(ia[0], ia[1] + 1) - sdf(ia[0], ia[1] - 1);
		real l = 0, m = 0, r = 0;
		if (abs(DX) > abs(DY)) {
			for (int i = -S; i < S; i++) {
				l += length_fraction(sdf(ia[0] + i, ia[1] - 1), sdf(ia[0] + i + 1, ia[1] - 1));
				m += length_fraction(sdf(ia[0] + i, ia[1]), sdf(ia[0] + i + 1, ia[1]));
				r += length_fraction(sdf(ia[0] + i, ia[1] + 1), sdf(ia[0] + i + 1, ia[1] + 1));
			}
		} else {
			for (int i = -S; i < S; i++) {
				l += length_fraction(sdf(ia[0] - 1, ia[1] + i), sdf(ia[0] - 1, ia[1] + i + 1));
				m += length_fraction(sdf(ia[0], ia[1] + i), sdf(ia[0], ia[1] + i + 1));
				r += length_fraction(sdf(ia[0] + 1, ia[1] + i), sdf(ia[0] + 1, ia[1] + i + 1));
			}
		}
		real curv = -hf_curv2d(l, m, r);
		ret[ia] = std::clamp<real>(curv, -2.0, 2.0);
	}
	return ret;
}

// return the curvature of a 3D level set, for cells near the surface.
template<typename real>
inline Array<3, real> getCurvature(const Array<3, real> &sdf, Array<3, int> surface) {
	using namespace range_operator;
	Array<3, real> ret(sdf.dim);

	ret = 0;
	int S = 3;
	const int T = 1;
	for (auto ia : Subrange<3>(sdf.dim, { { S, -S },{ S, -S },{ S, -S } })) {
		if (!surface[ia])
			continue;

		real DX = sdf(ia[0] + 1, ia[1], ia[2]) - sdf(ia[0] - 1, ia[1], ia[2]);
		real DY = sdf(ia[0], ia[1] + 1, ia[2]) - sdf(ia[0], ia[1] - 1, ia[2]);
		real DZ = sdf(ia[0], ia[1], ia[2] + 1) - sdf(ia[0], ia[1], ia[2] - 1);
		real L = sqrt(DX * DX + DY * DY + DZ * DZ);
		DX /= L; DY /= L; DZ /= L;
		real max;
		int orientation;
		if (abs(DX) > abs(DY)) {
			max = DX;
			orientation = 0;
		} else {
			max = DY;
			orientation = 1;
		}
		orientation = (abs(max) > abs(DZ)) ? orientation : 2;
		real theta = acos(abs(max));
		real gamma = theta < 0.8 ? 0.0 : 0.2;
		real hs[2 * T + 1][2 * T + 1];
		if (orientation == 0) {
			for (int i = -T; i < T + 1; i++) {
				for (int j = -T; j < T + 1; j++) {
					hs[i + T][j + T] = 0;
					for (int k = -S; k < S; k++) {
						hs[i + T][j + T] += length_fraction(
							sdf(ia[0] + k, ia[1] + i, ia[2] + j),
							sdf(ia[0] + k + 1, ia[1] + i, ia[2] + j)
						);
					}
				}
			}
		} else if (orientation == 1) {
			for (int i = -T; i < T + 1; i++) {
				for (int j = -T; j < T + 1; j++) {
					hs[i + T][j + T] = 0;
					for (int k = -S; k < S; k++) {
						hs[i + T][j + T] += length_fraction(
							sdf(ia[0] + i, ia[1] + k, ia[2] + j),
							sdf(ia[0] + i, ia[1] + k + 1, ia[2] + j)
						);
					}
				}
			}
		} else if (orientation == 2) {
			for (int i = -T; i < T + 1; i++) {
				for (int j = -T; j < T + 1; j++) {
					hs[i + T][j + T] = 0;
					for (int k = -S; k < S; k++) {
						hs[i + T][j + T] += length_fraction(
							sdf(ia[0] + i, ia[1] + j, ia[2] + k),
							sdf(ia[0] + i, ia[1] + j, ia[2] + k + 1)
						);
					}
				}
			}
		} else {
			std::cout << "ERR" << std::endl;
			exit(0);
		}
		real curv = -hf_curv3d(hs, gamma);
		ret[ia] = std::clamp<real>(curv, -2.0, 2.0);
	}
	return ret;
}

// changes value of @sdf for cells away from @surface, up to maximum band size, using fast march method
// @sdf values for cells farther than @narrow_band are truncated to (+-) @narrow_band
template<int N, typename real>
inline void fastMarching(Array<N, real> &sdf, real narrow_band, const Array<N, int> &surface) {
	using namespace range_operator;
	const int UNKNOWN = 0, INQUEUE = 1, KNOWN = 2;
	Array<N, int> tag = surface * 2;
	Array<N, real> usdf = sdf.map<real>([](real r) {return abs(r); });

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

	auto updater = [&usdf, &tag, &UNKNOWN, &INQUEUE, &KNOWN](auto &ia) -> real {
		real min = 65535, mins[N];
		real A = 0, B = 0, C = -1;
		for (int i = 0; i < N; i++) {
			mins[i] = 65535;
			if (ia[i] > 0) {
				if (tag[I<N>(ia, i, -1)] == KNOWN) {
					mins[i] = usdf[I<N>(ia, i, -1)];
				}
			}
			if (ia[i] + 1 < usdf.dim[i]) {
				if (tag[I<N>(ia, i, 1)] == KNOWN) {
					mins[i] = std::min<real>(mins[i], usdf[I<N>(ia, i, 1)]);
				}
			}
			if (mins[i] < 65534) {
				A += 1;
				B -= 2 * mins[i];
				C += mins[i] * mins[i];
			}
			min = std::min<real>(min, mins[i]);
		}
		real r = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
		r = (r < min + 1) ? r : min + 1;
		return r;
	};

	std::set<int, decltype(comparator)> Q(comparator);

	real mins[N];

	// initialize a priority queue
	auto offset = dim2offset<N>(sdf.dim);
	for (int j = 0; j < sdf.size(); j++) {
		if (tag[j] == KNOWN) {
			auto ia = i2ia<N>(j, offset);
			auto neighbors = usdf.getNeighborsOf_O2(ia);
			for (int k = 0; true; k++) {
				int nb = neighbors[k];
				if (nb < 0) break;
				int tag_nb = tag[nb];
				if (tag_nb != KNOWN) {
					tag[nb] = INQUEUE;
				}
			}
		} else {
			// is this necessary? todo
			sdf[j] = sdf[j] > 0 ? narrow_band : -narrow_band;
		}
	}

	for (int j = 0; j < sdf.size(); j++) {
		if (tag[j] == INQUEUE) {
			// compute new usdf
			real r = updater(i2ia<N>(j, offset));
			usdf[j] = r;
			Q.insert(j);
		}
	}

	while (!Q.empty()) {
		int p = *Q.begin();
		sdf[p] = sdf[p] > 0 ? usdf[p] : -usdf[p];
		tag[p] = KNOWN;
		if (usdf[p] > narrow_band)
			break;
		Q.erase(p);

		int ia[N];
		tag._getIndex(p, ia);

		auto neighbors = usdf.getNeighborsOf_O2(ia);
		for (int k = 0; true; k++) {
			int nb = neighbors[k];

			if (nb < 0)
				break;

			if (tag[nb] == KNOWN)
				continue;

			int nb_ia[N];
			tag._getIndex(nb, nb_ia);

			real r = updater(nb_ia);

			if (tag[nb] == UNKNOWN) {
				usdf[nb_ia] = r;
				tag[nb_ia] = INQUEUE;
				Q.insert(nb);
			} else if (tag[nb] == INQUEUE) {
				// pop and re enqueue
				if (r < usdf[nb]) {
					Q.erase(nb);
					usdf[nb] = r;
					Q.insert(nb);
				}
			}
		}
	}

	for (auto ia : Range<N>(sdf.dim)) {
		if (tag[ia] != KNOWN) {
			sdf[ia] = sdf[ia] > 0 ? narrow_band : -narrow_band;
		}
	}
}

// changes @sdf where @surface is true.
// new values are shortest distance to zero-interpolant along axes.
// returns those zero-interpolant positions
template<int N, typename real>
std::vector<std::array<real, N>> initialize_a(Array<N, real> &sdf, const Array<N, int> &surface) {
	using namespace range_operator;

	auto offset = dim2offset<N>(sdf.dim);
	std::vector<std::array<real, N>> sp;
	std::vector<real> sdf_update;

	for (int i = 0; i < sdf.size(); i++) {
		if (!surface[i])
			continue;
		auto ia = i2ia<N>(i, offset);
		real axis = -1;
		real d = 2;
		for (int j = 0; j < N; j++) {
			if (ia[j] > 0) {
				real s_l = sdf[I<N>(ia, j, -1)];
				if (!same_sign(s_l, sdf[i])) {
					real d_l = sdf[i] / (sdf[i] - s_l);
					if (d_l < abs(d)) {
						axis = j;
						d = -d_l;
					}
				}
			}
			if (ia[j] + 1 < sdf.dim[j]) {
				real s_r = sdf[I<N>(ia, j, 1)];
				if (!same_sign(s_r, sdf[i])) {
					real d_r = sdf[i] / (sdf[i] - s_r);
					if (d_r < abs(d)) {
						axis = j;
						d = d_r;
					}
				}
			}

		}
		std::array<real, N> p;
		copy_indexable<N>(ia, p);
		p[axis] += d;
		sp.push_back(p);
		sdf_update.push_back(sdf[i] >= 0 ? abs(d) : -abs(d));
	}
	int cursor = 0;
	for (int i = 0; i < sdf.size(); i++) {
		if (surface[i]) {
			sdf[i] = sdf_update[cursor++];
		}
	}
	return sp;
}

// return the intersection points of 0-contour and the grid
// no side effects on inputs
template<int N, typename real>
std::vector<std::array<real, N>> initialize_b(const Array<N, real> &sdf, const Array<N, int> &surface) {
	using namespace range_operator;
	std::vector<std::array<real, N>> sp;

	auto offset = dim2offset<N>(sdf.dim);
	for (int i = 0; i < sdf.size(); i++) {
		if (!surface[i])
			continue;

		auto ia = i2ia<N>(i, offset);

		if (sdf[i] == 0) {
			std::array<real, N> p;
			copy_indexable<N>(ia, p);
			sp.push_back(p);
			continue;
		}

		for (int i = 0; i < N; i++) {
			if (ia[i] + 1 == sdf.dim[i])
				continue;
			auto ia_r = I<N>(ia, i, 1);
			if (sdf[ia] * sdf[ia_r] < 0) {
				std::array<real, N> p;
				copy_indexable<N>(ia, p);
				real theta = sdf[ia] / (sdf[ia] - sdf[ia_r]);
				p[i] += theta;
				sp.push_back(p);
			}
		}
	}

	return sp;
}

// @sp is a set of surface points returned from either initialize_a() or initialize_b()
// @sdf within @band distance from @sp are updated, as their shortest distance to the point set @sp
// returns: @i2sp[@ia] = -1 if @ia is outside the band, otherwise @sp[@i2sp[@ia]] is the point in @sp that's cloesest to @ia.
template<int N, typename real>
Array<N, int> gi2sp(Array<N, real> &sdf, const std::vector<std::array<real, N>> &sp, real band) {
	using namespace range_operator;
	real band_sq = band * band;
	Array<N, int> i2sp(sdf.dim, -1);
	Array<N, real> dsq(sdf.dim, band_sq);
	auto D_SQ = [](auto x, auto y) {
		real sum = 0;
		for (int i = 0; i < N; i++) {
			sum += (x[i] - y[i]) * (x[i] - y[i]);
		}
		return sum;
	};

	for (int i = 0; i < sp.size(); i++) {
		int low[N], r[N];
		for (int j = 0; j < N; j++) {
			low[j] = std::max<int>((int)(sp[i][j] - band) + 1, 0);
			int hi = std::min<int>(sp[i][j] + band, sdf.dim[j]);
			r[j] = hi - low[j];
		}
		for (auto r_ia : Range<N>(r)) {
			int ia[N];
			for (int j = 0; j < N; j++) {
				ia[j] = r_ia[j] + low[j];
			}
			real dsq1 = D_SQ(ia, sp[i]);
			if (dsq1 < dsq[ia]) {
				dsq[ia] = dsq1;
				i2sp[ia] = i;
			}
		}
	}

	for (int i = 0; i < sdf.size(); i++) {
		if (i2sp[i] < 0) {
			sdf[i] = sdf[i] >= 0 ? band : -band;
		} else {
			sdf[i] = sdf[i] >= 0 ? sqrt(dsq[i]) : -sqrt(dsq[i]);
		}
	}

	return i2sp;
}
