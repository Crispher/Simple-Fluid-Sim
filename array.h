#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cassert>
#include "range_operator.h"

// todo:
// range for should replace forEach()
// get_functor should be replaced
// add begin() end() to Array/SubArray

/* type traits*/
template<int N, typename T>
class Array;

template<int N, typename T>
class SubArray;

template <typename T>
struct is_indexable_ints : std::conditional<
	std::is_convertible_v<T, const int* const> || 
	std::is_convertible_v<T, const std::vector<int>> 
	, 
	std::true_type,
	std::false_type
>::type {};


template <typename T, int N>
struct is_array_ints : std::conditional<
	std::is_convertible_v<T, const std::array<int, N>>
	,
	std::true_type,
	std::false_type
>::type {
};


/***********************
	CLASS SUBARRAY
***********************/

template<int N, typename T>
class SubArray {
	template <int N, typename T2>
	friend class Array;
private:
	Array<N, T>* pArray;
	int lo[N];
	int hi[N];
	int dim[N];
	bool is_degenerate[N]; // todo

private:
	SubArray(Array<N, T> *p, int *lo, int *hi, int *dim, bool *degenerate);
public:
	~SubArray() {}

public:
	void operator=(const T&t) { forEach([&t](int*, T &v) { v = t; }); }
	void operator+=(const T&t) { forEach([&t](int*, T &v) { v += t; }); }
	void operator-=(const T&t) { forEach([&t](int*, T &v) { v -= t; }); }
	void operator*=(const T&t) { forEach([&t](int*, T &v) { v *= t; }); }
	void operator/=(const T&t) { forEach([&t](int*, T &v) { v /= t; }); }
	
	void operator=(const Array<N, T> &t) { forEach([&t](int* i, T&v) { v = t[i]; }); }
	void operator=(const SubArray<N, T> &t) { operator=(Array<N, T>(t)); }
	void operator+=(const Array<N, T>&t) { forEach([&t](int* i, T&v) { v += t[i]; }); }
	void operator-=(const Array<N, T>&t) { forEach([&t](int* i, T&v) { v -= t[i]; }); }
	void operator*=(const Array<N, T>&t) { forEach([&t](int* i, T&v) { v *= t[i]; }); }
	void operator/=(const Array<N, T>&t) { forEach([&t](int* i, T&v) { v /= t[i]; }); }

	Array<N, T> operator+(const T &t) const { return Array<N, T>(*this) + t; }
	Array<N, T> operator-(const T &t) const { return Array<N, T>(*this) - t; }
	Array<N, T> operator*(const T &t) const { return Array<N, T>(*this) * t; }
	Array<N, T> operator/(const T &t) const { return Array<N, T>(*this) / t; }

	Array<N, T> operator+(const Array<N, T> &t) const { return Array<N, T>(*this) + t; }
	Array<N, T> operator-(const Array<N, T> &t) const { return Array<N, T>(*this) - t; }
	Array<N, T> operator*(const Array<N, T> &t) const { return Array<N, T>(*this) * t; }
	Array<N, T> operator/(const Array<N, T> &t) const { return Array<N, T>(*this) / t; }
	
public:
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value>>
	T operator[] (Indexable ia) const { return (*pArray)[_getLocation(ia)]; }
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value>> 
	T& operator[] (Indexable ia) { return (*pArray)[_getLocation(ia)]; }

private:
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value>>
	int _getLocation(Indexable) const;
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value>>
	T _get(Indexable ia) const;
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value>>
	T& _get(Indexable ia);
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value>>
	bool _next(Indexable &ia) const;

	void forEach(std::function<void(int*, T)>) const;
	void forEach(std::function<void(int*, T&)>);


};

template<int N, typename T>
SubArray<N, T>::SubArray(Array<N, T> *p, int *lo, int *hi, int *dim, bool *degenerate) {
	pArray = p;
	for (int i = 0; i < N; i++) {
		this->lo[i] = lo[i];
		this->hi[i] = hi[i];
		this->dim[i] = dim[i];
		this->is_degenerate[i] = degenerate[i];
	}
}

template<int N, typename T>
template<typename Indexable, typename check>
inline int SubArray<N, T>::_getLocation(Indexable ia) const {
	int ind[N];
	for (int i = 0; i < N; i++) {
		ind[i] = lo[i] + ia[i];
	}
	return pArray->_getLocation(ind);
}

template<int N, typename T>
template<typename Indexable, typename check>
inline T SubArray<N, T>::_get(Indexable ia) const {
	int ind[N];
	for (int i = 0; i < N; i++) {
		ind[i] = lo[i] + ia[i];
	}
	return (*pArray)[ind];
}

template<int N, typename T>
template<typename Indexable, typename check>
inline T & SubArray<N, T>::_get(Indexable ia) {
	static_assert(!std::is_const_v<decltype(*pArray)>);
	int ind[N];
	for (int i = 0; i < N; i++) {
		ind[i] = lo[i] + ia[i];
	}
	return (*pArray)[ind];
}


template<int N, typename T>
template<typename Indexable, typename check>
inline bool SubArray<N, T>::_next(Indexable &curr) const {
	for (int i = N - 1; i > 0; i--) {
		curr[i]++;
		if (curr[i] == dim[i]) {
			curr[i] = 0;
		} else {
			return true;
		}
	}
	curr[0]++;
	return curr[0] < dim[0];
}

template<int N, typename T>
inline void SubArray<N, T>::forEach(std::function<void(int*, T)> func) const {
	int curr[N] = { 0 };
	do {
		func(curr, _get(curr));
	} while (_next(curr));
}

template<int N, typename T>
inline void SubArray<N, T>::forEach(std::function<void(int*, T&)> func) {
	int curr[N] = { 0 };
	do {
		func(curr, _get(curr));
	} while (_next(curr));
}


/*************************
	ARRAY CLASS TEMPLATE
*************************/


// to be replaced by Operator interface
template<int N, typename R, typename T1, typename ... Ts>
std::function<Array<N, R>(const Array<N, T1>&, const Array<N, Ts>& ...)> get_array_functor(std::function<R(T1, Ts...)> f) {
	return [f](const Array<N, T1> &a1, const Array<N, Ts>& ... as) -> Array<N, R> {
		return Array<N, R>([&](int *ia) {return f(a1[ia], as[ia]...); }, a1.dim);
	};
}

namespace array_util {
	typedef struct {} CopyOutermost;
	inline static CopyOutermost *COPY_OUTERMOST = nullptr;

}

template<int N, typename T>
class Array {
	template<int N, typename T2>
	friend class Array;
public: // accessor
	const int *dim = _dim;
	int size() const { return _size; }

private:
	int _dim[N];
	int _size;
	std::vector<T> data;
	int offset[N];

private:
	bool _initialized = true;

public:
	Array() { _initialized = false; }

	Array(std::ifstream& is) { load(is); }

	// 4 constructors, tested
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value>>
	Array(Indexable dims);
	Array(std::initializer_list<int>);
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value>>
	Array(Indexable dims, const T& t);

	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value >>
	Array(std::function<T(const int*)>, Indexable);
	Array(std::function<T(const int*)>, std::initializer_list<int>);

	template<typename T2, typename check = std::enable_if_t<std::is_convertible_v<T2, T>>>
	Array(const SubArray<N, T2>&);

	Array(const Array<N, T>& t);

	template<int N2, typename check = std::enable_if_t<(N2 > N)>>
	Array(const SubArray<N2, T> &t);

	~Array() {};

private:
	Array<N, T>& _copy(const Array<N, T> &t) {
		if (!t._initialized) {
			printf("copy of uninitialized array\n");
			_initialized = false;
			return *this;
		}
		for (int i = 0; i < N; i++) {
			_dim[i] = t._dim[i];
			offset[i] = t.offset[i];
		}
		data = t.data;
		_size = t._size;
		_initialized = true;
		return *this;
	}
public:

//private:
	// the only interface to compute location from index
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value>>
	int _getLocation(Indexable) const;

	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value>>
	bool _getIndex(int loc, Indexable&) const;

	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value>>
	T& _get(Indexable);

	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value>>
	T _get(Indexable) const;

	// increment the index, return true iff the incremented index lies in boundary
	// Indexable can be int* or std::vector<int>
	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value >>
	bool _next(Indexable&) const;

	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value || is_array_ints<Indexable, N>::value >>
	bool _checkIndex(const Indexable&) const;

	void _c(const T& = T()) const {
		assert(_initialized && "Use of uninitialized array.");
	}

	template<typename T2>
	void _c(const Array<N, T2>& t) const {
		assert(_initialized && t._initialized && "Use of uninitialized array.");
		assert(check_dim_agree(*this, t) && "Dimension does not agree.");
	}

public:
	// indexing operator, taking 1 int, return data[loc]
	T operator[](int loc) const { return data[loc]; }
	T& operator[](int loc) { return data[loc]; }

	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value||is_array_ints<Indexable, N>::value>>
	T operator[](Indexable) const;
	T operator[](std::initializer_list<int>) const;

	template<typename Indexable, typename check = std::enable_if_t<is_indexable_ints<Indexable>::value||is_array_ints<Indexable, N>::value>>
	T& operator[](Indexable);
	T& operator[](std::initializer_list<int>);

	// indexing operator, taking N_DIM ints, manual interface
	template<typename ... types, typename check=int>// = std::enable_if_t<sizeof...(types) == N_DIM - 1, int>>
	T& operator()(int i, types ... args);
	template<typename ... types, typename check=int>// = std::enable_if_t<sizeof...(types) == N_DIM - 1, int>>
	T operator()(int i, types ... args) const;

public:
	// slicing operator
	const SubArray<N, T> operator[](int r[N][2]) const;
	const SubArray<N, T> operator[](std::vector<std::vector<int>>) const;
	SubArray<N, T> operator[](int r[N][2]);
	SubArray<N, T> operator[](std::vector<std::vector<int>>);

public:

	// operator defined for all array of T2 convertible to T, and subarray of type T2 convertible to T, via implicit conversion
	
	Array<N, T>& operator=(const T& t) { _c(t); for (T&d : data) d = t; return (*this); }
	Array<N, T>& operator+=(const T& t) { _c(t); for (T&d : data) d += t; return (*this); }
	Array<N, T>& operator-=(const T& t) { _c(t); for (T&d : data) d -= t; return (*this); }
	Array<N, T>& operator*=(const T& t) { _c(t); for (T&d : data) d *= t; return (*this); }
	Array<N, T>& operator/=(const T& t) { _c(t); for (T&d : data) d /= t; return (*this); }
	
	Array<N, T>& operator=(const Array<N, T>& t) { if (!_initialized) return _copy(t);
														   _c(t); for (int i = 0; i < _size; i++) data[i] = t[i]; return(*this); }
	Array<N, T>& operator+=(const Array<N, T>&t) { _c(t); for (int i = 0; i < _size; i++) data[i] += t[i]; return(*this); }
	Array<N, T>& operator-=(const Array<N, T>&t) { _c(t); for (int i = 0; i < _size; i++) data[i] -= t[i]; return(*this); }
	Array<N, T>& operator*=(const Array<N, T>&t) { _c(t); for (int i = 0; i < _size; i++) data[i] *= t[i]; return(*this); }
	Array<N, T>& operator/=(const Array<N, T>&t) { _c(t); for (int i = 0; i < _size; i++) data[i] /= t[i]; return(*this); }

	Array<N, T> operator+(const T& t) const { _c(t); return Array<N, T>(*this) += t; }
	Array<N, T> operator-(const T& t) const { _c(t); return Array<N, T>(*this) -= t; }
	Array<N, T> operator*(const T& t) const { _c(t); return Array<N, T>(*this) *= t; }
	Array<N, T> operator/(const T& t) const { _c(t); return Array<N, T>(*this) /= t; }

	Array<N, T> operator+(const Array<N, T>& t) const { _c(t); return Array<N, T>(*this) += t; }
	Array<N, T> operator-(const Array<N, T>& t) const { _c(t); return Array<N, T>(*this) -= t; }
	Array<N, T> operator*(const Array<N, T>& t) const { _c(t); return Array<N, T>(*this) *= t; }
	Array<N, T> operator/(const Array<N, T>& t) const { _c(t); return Array<N, T>(*this) /= t; }

	typedef std::function<T(const T&, const T&)> _bin_op_t;
	Array<N, T> cwise(_bin_op_t op, const T& rhs) 				const { Array<N, T> ret(dim); for (int i = 0; i < _size; i++) ret[i] = op(data[i], rhs); return ret; }
	Array<N, T> cwise(_bin_op_t op, const Array<N, T> &rhs) const { Array<N, T> ret(dim); for (int i = 0; i < _size; i++) ret[i] = op(data[i], rhs[i]); return ret; }


public:
	// implicit conversion
	template<typename T2, typename check = std::enable_if_t<std::is_convertible_v<T, T2>>>
	operator Array<N, T2>() const { 
		Array<N, T2> ret(_dim);
		for (int i = 0; i < _size; i++) {
			ret.data[i] = data[i];
		}
		return ret;
	}

public:
	// in the following 2 functions: a non-grid point P is clamped s.t. along all axes,
	// P[i] \in [0, dim[i] - 1)

	// return interpolated value at a non-grid point
	template<typename IndexableReal>
	T interpolate(const IndexableReal &) const;
	// return weighted interpolated value at a non-grid point, RAISES ERROR if sum(weight) == 0;
	template<typename IndexableReal, typename WeightFactorType>
	T interpolate(const IndexableReal &, const Array<N, WeightFactorType>&) const;
	template<typename IndexableReal, bool clamp = true>
	T interpolate_cubic(const IndexableReal &) const;

public:
	void forEach(std::function<void(int*, T&)>);
	void forEach(std::function<void(int*, T)>) const;
	std::vector<T> flatten() const;
	void setData(const std::vector<T>&d) { data = d; }
	template<typename Indexable>
	int getNeighbor(const Indexable& ia, int axis, int distance) const;
	template<typename Indexable>
	std::vector<int> getNeighborsOf(const Indexable&) const;
	template<typename Indexable>
	std::array<int, 2*N + 1> getNeighborsOf_O2(const Indexable&) const;
	template<typename Indexable>
	std::vector<int> getNeighborsOf(const Indexable&, int along_axis) const;

public:
	template<typename Target>
	Array<N, Target> map(std::function<Target(T)>) const;
	Array<N, T> pad(std::vector<std::vector<int>>, const T&v) const;
	Array<N, T> pad(std::vector<std::vector<int>>, array_util::CopyOutermost *c) const;
	Array<N, T> pad(int axis, int l, int r, const T&v) const;
	Array<N, T> pad(int axis, int l, int r, array_util::CopyOutermost *c) const;

public:
	void print(std::string name = "") const;
	void save(std::ostream &os) const;
	void load(std::istream &is);
	T max() const { return *std::max_element(data.begin(), data.end()); }
	T max(std::function<T(T)> f) { T M = f(data[0]); for (T &t : data) M = std::max<T>(f(t), M); return M; }
	T min() const { return *std::min_element(data.begin(), data.end()); }
	T sum() const { T ret = 0; for (auto &i : data) { ret += i; } return ret; }
	bool isInitialized() const { return _initialized; }
	template<typename Indexable>
	bool isInRange(const Indexable &ia) { return _checkIndex(ia); }
	void disableWarning() { _initialized = true; }
	bool hasNaN() const { return std::any_of<decltype(data.begin()), bool(T)>(data.begin(), data.end(), std::isnan<T>); }
};


#pragma region CONSTRUCTORS

template<int N, typename T>
template<typename Indexable, typename check>
inline Array<N, T>::Array(Indexable dims) {
	_size = 1;
	for (int i = 0; i < N; i++) {
		_dim[i] = dims[i];
		_size *= _dim[i];
	}
	offset[N - 1] = 1;
	for (int i = N - 2; i >= 0; i--) {
		offset[i] = _dim[i + 1] * offset[i + 1];
	}
	data = std::vector<T>(_size);
}

template<int N, typename T>
template<typename Indexable, typename check>
inline Array<N, T>::Array(Indexable dims, const T & t) : Array<N, T>(dims) {
	for (auto &e : data) {
		e = t;
	}
}

template<int N, typename T>
inline Array<N, T>::Array(std::initializer_list<int> dims) :
	Array(std::vector<int>(dims)) {
}

template<int N, typename T>
template<typename Indexable, typename check>
inline Array<N, T>::Array(std::function<T(const int*)> func, Indexable dims) :
	Array(dims) {
	using namespace range_operator;
	for (auto &[ia, val] : enumerate(*this)) {
		val = func(ia);
	}
}

template<int N, typename T>
inline Array<N, T>::Array(std::function<T(const int*)> func, std::initializer_list<int> dims) : 
	Array<N, T>(func, std::vector<int>(dims)) {
}

template<int N, typename T>
template<typename T2, typename check>
inline Array<N, T>::Array(const SubArray<N, T2>& s) :
	Array(s.dim) {
	forEach([&s](int *p, T& val) {
		val = s[p];
	});
}

template<int N, typename T>
inline Array<N, T>::Array(const Array<N, T>& t) {
	for (int i = 0; i < N; i++) {
		_dim[i] = t._dim[i];
		offset[i] = t.offset[i];
	}
	data = t.data;
	_size = t._size;
	_initialized = true;
}


template<int N, typename T>
template<int N2, typename check>
inline Array <N, T>::Array(const SubArray<N2, T>& t) {
	int j = 0;
	_size = 1;
	for (int i = 0; i < N2; i++) {
		if (!t.is_degenerate[i]) {
			_dim[j++] = t.dim[i];
			_size *= t.dim[i];
		}
	}
	assert(j == N && "cannot convert subarray to lower dimension");
	offset[N - 1] = 1;
	for (int i = N - 2; i >= 0; i--) {
		offset[i] = _dim[i + 1] * offset[i + 1];
	}
	data = (Array<N2, T>(t)).data;
}

#pragma endregion

#pragma region HELPERS
template<int N, typename T>
template<typename Indexable, typename check>
inline int Array<N, T>::_getLocation(Indexable ia) const {
	int loc = 0;
	for (int i = 0; i < N; i++) {
		loc += ia[i] * offset[i];
	}
	return loc;
}

template<int N, typename T>
template<typename Indexable, typename check>
inline bool Array<N, T>::_getIndex(int loc, Indexable &ia) const {
	ia[0] = loc / offset[0];
	for (int i = 1; i < N; i++) {
		ia[i] = (loc % offset[i - 1]) / offset[i];
	}
	return (0 <= loc && loc < _size);
}

template<int N, typename T>
template<typename Indexable, typename check>
inline T & Array<N, T>::_get(Indexable ia) {
	return data[_getLocation(ia)];
}

template<int N, typename T>
template<typename Indexable, typename check>
inline T Array<N, T>::_get(Indexable ia) const{
	return data[_getLocation(ia)];
}

template<int N, typename T>
template<typename Indexable, typename check>
inline bool Array<N, T>::_next(Indexable &curr) const {
	for (int i = N - 1; i > 0; i--) {
		curr[i]++;
		if (curr[i] == _dim[i]) {
			curr[i] = 0;
		} else {
			return true;
		}
	}
	curr[0]++;
	return curr[0] < _dim[0];
}

template<int N, typename T>
template<typename Indexable, typename check>
inline bool Array<N, T>::_checkIndex (const Indexable &ia) const {
	for (int i = 0; i < N; i++) {
		if (ia[i] < 0 || ia[i] >= _dim[i]) {
			return false;
		}
	}
	return true;
}

#pragma endregion

template<int N, typename T>
template<typename Indexable, typename check>
inline T Array<N, T>::operator[](Indexable a) const {
	return data[_getLocation(a)];
}


template<int N, typename T>
inline T Array<N, T>::operator[](std::initializer_list<int> a) const {
	return data[_getLocation(std::vector<int>(a))];
}

template<int N, typename T>
template<typename Indexable, typename check>
inline T& Array<N, T>::operator[](Indexable a) {
	return _get(a);
}


template<int N, typename T>
inline T& Array<N, T>::operator[](std::initializer_list<int> a) {
	return _get(std::vector<int>(a));
}

template<int N, typename T>
template<typename ...types, typename check>
inline T Array<N, T>::operator()(int i, types ...args) const {
	int ia[N] = { i, args... };
	return _get(ia);
}

template<int N, typename T>
template<typename ...types, typename check>
inline T& Array<N, T>::operator()(int i, types ...args) {
	int ia[N] = { i, args... };
	return _get(ia);
}

template<int N, typename T>
template<typename IndexableReal>
inline T Array<N, T>::interpolate(const IndexableReal &v) const {
	using namespace range_operator;
	int lo[N], hi[N];
	double clamped[N];
	int cellIndex[N][2];
	double wa[N][2];
	int twos[N];
	for (int i = 0; i < N; i++) {
		twos[i] = 2;
	}

	for (int i = 0; i < N; i++) {
		clamped[i] = std::clamp<double>(v[i], 0.0, _dim[i] - 1.0);
		lo[i] = std::min<int>((int)(clamped[i]), _dim[i] - 2);
		hi[i] = lo[i] + 1;
		cellIndex[i][0] = lo[i]; 
		cellIndex[i][1] = lo[i] + 2;
		wa[i][0] = hi[i] - clamped[i];
		wa[i][1] = 1 - wa[i][0];
	}

	std::array<double, 1 << N> p_weight;
	std::array<double, 1 << N> cell_a;
	for (auto ia : Range<N>(twos)) {
		p_weight[ia[N]] = 1.0;
		for (int i = 0; i < N; i++) {
			p_weight[ia[N]] *= wa[i][ia[i]];
		}

		int id[N];
		for (int i = 0; i < N; i++) {
			id[i] = ia[i] + cellIndex[i][0];
		}
		cell_a[ia[N]] = (*this)[id];
	}
	
	T ret2 = 0;
	for (int i = 0; i < (1 << N); i++) {
		ret2 += cell_a[i] * p_weight[i];
	}

	return ret2;
}

template<int N, typename T>
template<typename IndexableReal, typename WeightFactorType>
inline T Array<N, T>::interpolate(const IndexableReal &v, const Array<N, WeightFactorType>& wf) const {
	using namespace range_operator;
	int lo[N], hi[N];
	double clamped[N];
	std::vector<std::vector<int>> cellIndex(N);
	double wa[N][2];
	int twos[N];
	for (int i = 0; i < N; i++) {
		twos[i] = 2;
	}
	for (int i = 0; i < N; i++) {
		clamped[i] = std::clamp<double>(v[i], 0, _dim[i] - 1.0);
		lo[i] = std::min<int>((int)(clamped[i]), _dim[i] - 2);
		hi[i] = lo[i] + 1;
		cellIndex[i] = { lo[i], lo[i] + 2 };
		wa[i][0] = hi[i] - clamped[i];
		wa[i][1] = 1 - wa[i][0];
	}


	std::array<double, 1 << N> p_weight;
	for (auto ia : Range<N>(twos)) {
		p_weight[ia[N]] = 1.0;
		for (int i = 0; i < N; i++) {
			p_weight[ia[N]] *= wa[i][ia[i]];
		}
		// p_weight[ia[N_DIM]] += 1e-5;
	}

	std::array<double, 1 << N> cell_a;
	std::array<double, 1 << N> wfcell_a;
	for (auto ia : Range<N>(twos)) {
		int id[N];
		for (int i = 0; i < N; i++) {
			id[i] = ia[i] + cellIndex[i][0];
		}
		cell_a[ia[N]] = (*this)[id];
		wfcell_a[ia[N]] = wf[id];
	}


	std::array<double, 1 << N> finalWeight_a;
	T weightedsum = 0, weightsum = 0;
	for (int i = 0; i < (1 << N); i++) {
		finalWeight_a[i] = p_weight[i] * wfcell_a[i];
		weightsum += finalWeight_a[i];
		weightedsum += cell_a[i] * finalWeight_a[i];
	}

	weightedsum /= weightsum;
	if (weightsum <= 0) {
		assert(false);
	}
	return weightedsum;
}

template<int N, typename T>
template<typename IndexableReal, bool clamp>
inline T Array<N, T>::interpolate_cubic(const IndexableReal &v) const {
	using namespace range_operator;

	for (int i = 0; i < N; i++) {
		if (v[i] <= 1 || dim[i] - 2 <= v[i]) {
			return interpolate(v);
		}
	}

	int lo[N], hi[N];
	double wa[N][4];
	int twos[N];
	int fours[N];
	for (int i = 0; i < N; i++) {
		twos[i] = 2;
		fours[i] = 4;
	}

	for (int i = 0; i < N; i++) {
		lo[i] = (int)(v[i]) - 1;
		hi[i] = lo[i] + 4;
		double p = v[i] - lo[i];
		wa[i][0] = (p - 1) * (p - 2) * (p - 3) / (0. - 1) / (0. - 2) / (0. - 3);
		wa[i][1] = (p - 0) * (p - 2) * (p - 3) / (1. - 0) / (1. - 2) / (1. - 3);
		wa[i][2] = (p - 0) * (p - 1) * (p - 3) / (2. - 0) / (2. - 1) / (2. - 3);
		wa[i][3] = (p - 0) * (p - 1) * (p - 2) / (3. - 0) / (3. - 1) / (3. - 2);
	}

	std::array<double, 1 << (2 * N)> p_weight;
	std::array<T, 1 << (2 * N)> cell;
	std::array<T, 1 << N> cell_c;

	for (auto ia : Range<N>(fours)) {
		p_weight[ia[N]] = 1.0;
		for (int i = 0; i < N; i++) {
			p_weight[ia[N]] *= wa[i][ia[i]];
		}
		int id[N];
		for (int i = 0; i < N; i++) {
			id[i] = ia[i] + lo[i];
		}
		cell[ia[N]] = (*this)[id];
	}
	for (auto ia : Range<N>(twos)) {
		int id[N];
		for (int i = 0; i < N; i++) {
			id[i] = ia[i] + lo[i] + 1;
		}
		cell_c[ia[N]] = (*this)[id];
	}

	T ret = 0, min = cell_c[0], max = cell_c[0];
	for (int i = 0; i < (1 << (2 * N)); i++) {
		ret += cell[i] * p_weight[i];
	}
	if (clamp) {
		for (int i = 0; i < (1 << N); i++) {
			min = std::min<T>(min, cell_c[i]);
			max = std::max<T>(max, cell_c[i]);
		}
		ret = std::clamp<T>(ret, min, max);
	}
	return ret;
}

template<int N, typename T>
template<typename Indexable>
inline int Array<N, T>::getNeighbor(const Indexable &ia, int axis, int distance) const {
	if (ia[axis] + distance < 0 || ia[axis] + distance >= dim[axis]) {
		return -1;
	}
	int loc = _getLocation(ia);
	loc += offset[axis] * distance;
	return loc;
}

template<int N, typename T>
template<typename Indexable>
inline std::vector<int> Array<N, T>::getNeighborsOf(const Indexable &a) const {
	int ia[N];
	std::vector<int> ret;
	for (int i = 0; i < N; i++) {
		ia[i] = a[i];
	}
	for (int i = 0; i < N; i++) {
		ia[i] --;
		if (_checkIndex(ia)) {
			ret.push_back(_getLocation(ia));
		}
		ia[i] += 2;
		if (_checkIndex(ia)) {
			ret.push_back(_getLocation(ia));
		}
		ia[i] --; // restore and go on to next axis;
	}
	return ret;
}

template<int N, typename T>
template<typename Indexable>
inline std::array<int, 2 * N + 1> Array<N, T>::getNeighborsOf_O2(const Indexable & a) const
{
	std::array<int, 2 * N + 1> ret;
	int cursor = 0;
	int ia[N];
	for (int i = 0; i < N; i++) {
		ia[i] = a[i];
	}
	for (int i = 0; i < N; i++) {
		ia[i] --;
		if (_checkIndex(ia)) {
			ret[cursor++] = _getLocation(ia);
		}
		ia[i] += 2;
		if (_checkIndex(ia)) {
			ret[cursor++] = _getLocation(ia);
		}
		ia[i] --; // restore and go on to next axis;
	}
	ret[cursor] = -1;
	return ret;
}

template<int N, typename T>
template<typename Indexable>
inline std::vector<int> Array<N, T>::getNeighborsOf(const Indexable &a, int along_axis) const {
	int ia[N];
	std::vector<int> ret;
	for (int i = 0; i < N; i++) {
		ia[i] = a[i];
	}
	
	const int i = along_axis;

	ia[i] --;
	if (_checkIndex(ia)) {
		ret.push_back(_getLocation(ia));
	}
	ia[i] += 2;
	if (_checkIndex(ia)) {
		ret.push_back(_getLocation(ia));
	}
	ia[i] --; // restore and go on to next axis;
	
	return ret;
}

template<int N, typename T>
template<typename Target>
inline Array<N, Target> Array<N, T>::map(std::function<Target(T)> f) const{
	Array<N, Target> ret(_dim);
	for (int i = 0; i < _size; i++) {
		ret.data[i] = f(data[i]);
	}
	return ret;
}

template<int N, typename T>
inline const SubArray<N, T> Array<N, T>::operator[](int r[N][2]) const {
	int lo[N], hi[N], sdim[N];
	bool degenerate[N];
	for (int i = 0; i < N; i++) {
		lo[i] = r[i][0];
		hi[i] = r[i][1];
		sdim[i] = hi[i] - lo[i];
		degenerate[i] = false;
	}
	return SubArray<N, T>(const_cast<Array<N, T>*>(this), lo, hi, sdim, degenerate);
}

template<int N, typename T>
inline const SubArray<N, T> Array<N, T>::operator[](std::vector<std::vector<int>> r) const {
	int lo[N], hi[N], sdim[N];
	bool degenerate[N];
	assert(r.size() == N);
	for (int i = 0; i < N; i++) {
		for (int p = 0; p < r[i].size(); p++) {
			if (r[i][p] < 0) {
				r[i][p] += _dim[i];
			}
			if (r[i][p] == underscore_as_semicolon::_) {
				r[i][p] = (p == r[i].size() - 1 ? _dim[i] : 0);
			}
		}
		if (r[i].size() == 0) {
			lo[i] = 0; hi[i] = _dim[i]; sdim[i] = _dim[i]; degenerate[i] = false;
		} else if (r[i].size() == 1) {
			lo[i] = r[i][0]; hi[i] = lo[i] + 1; sdim[i] = 1; degenerate[i] = true;
		} else if (r[i].size() == 2) {
			lo[i] = r[i][0]; hi[i] = r[i][1]; sdim[i] = hi[i] - lo[i]; degenerate[i] = false;
		} else {
			assert(false);
		}
	}
	return SubArray<N, T>(const_cast<Array<N, T>*>(this), lo, hi, sdim, degenerate);
	
}

template<int N, typename T>
inline SubArray<N, T> Array<N, T>::operator[](int r[N][2]) {
	int lo[N], hi[N], sdim[N];
	bool degenerate[N];
	for (int i = 0; i < N; i++) {
		lo[i] = r[i][0];
		hi[i] = r[i][1];
		sdim[i] = hi[i] - lo[i];
		degenerate[i] = false;
	}
	return SubArray<N, T>(this, lo, hi, sdim, degenerate);
}

template<int N, typename T>
inline SubArray<N, T> Array<N, T>::operator[](std::vector<std::vector<int>> r) {
	int lo[N], hi[N], sdim[N];
	bool degenerate[N];
	assert(r.size() == N);
	for (int i = 0; i < N; i++) {
		for (int p = 0; p < r[i].size(); p++) {
			if (r[i][p] < 0) {
				r[i][p] += _dim[i];
			}
			if (r[i][p] == underscore_as_semicolon::_) {
				r[i][p] = (p == r[i].size() - 1 ? _dim[i] : 0);
			}
		}
		if (r[i].size() == 0) {
			lo[i] = 0; hi[i] = _dim[i]; sdim[i] = _dim[i]; degenerate[i] = false;
		} else if (r[i].size() == 1) {
			lo[i] = r[i][0]; hi[i] = lo[i] + 1; sdim[i] = 1; degenerate[i] = true;
		} else if (r[i].size() == 2) {
			lo[i] = r[i][0]; hi[i] = r[i][1]; sdim[i] = hi[i] - lo[i]; degenerate[i] = false;
		} else {
			assert(false);
		}
	}
	return SubArray<N, T>(this, lo, hi, sdim, degenerate);
}

template<int N, typename T>
inline void Array<N, T>::forEach(std::function<void(int*, T&)> func) {
	int curr[N] = {0};
	do {
		func(curr, _get(curr));
	} while (_next(curr));
}

template<int N, typename T>
inline void Array<N, T>::forEach(std::function<void(int*, T)> func) const {
	int curr[N] = { 0 };
	do {
		func(curr, _get(curr));
	} while (_next(curr));
}

template<int N, typename T>
inline std::vector<T> Array<N, T>::flatten() const {
	return data;
}

template<int N, typename T>
inline Array<N, T> Array<N, T>::pad(std::vector<std::vector<int>> p, const T&v) const {
	int newdim[N];
	for (int i = 0; i < N; i++) {
		newdim[i] = p[i][0] + _dim[i] + p[i][1];
		p[i][1] = p[i][0] + _dim[i];
	}
	Array<N, T> ret(newdim);
	ret = v;
	ret[p] = (*this);
	return ret;
}

template<int N, typename T>
inline Array<N, T> Array<N, T>::pad(std::vector<std::vector<int>> p, array_util::CopyOutermost * c) const {
	Array<N, T> ret[N];
	ret[0] = this->pad(0, p[0][0], p[0][1], c);
	for (int i = 1; i < N; i++) {
		ret[i] = ret[i - 1].pad(i, p[i][0], p[i][1], c);
	}
	return ret[N-1];
}

template<int N, typename T>
inline Array<N, T> Array<N, T>::pad(int axis, int l, int r, const T & v) const {
	std::vector<std::vector<int>> padding(N, { 0, 0 });
	padding[axis] = { l, r };
	return pad(padding, v);
}

template<int N, typename T>
inline Array<N, T> Array<N, T>::pad(int axis, int l, int r, array_util::CopyOutermost *c) const {
	using namespace range_operator;
	Array<N, T> ret = pad(axis, l, r, T());
    std::vector<std::vector<int>> range = get_full_range<N>(), lb = range, rb = range;
	lb[axis] = { 0 }; rb[axis] = { -1 };
	for (int i = 0; i < l; i++) {
		range[axis] = { i };
		ret[R<N>(axis, { i })] = operator[](lb);
	}
	for (int i = l + _dim[axis]; i < l + _dim[axis] + r; i++) {
		range[axis] = { i };
		ret[R<N>(axis, { i })] = operator[](rb);
	}
	return ret;
}

template<int N, typename T>
inline void Array<N, T>::print(std::string name) const {
	if (name.size() > 0) {
		std::cout << name << " = \n";
	}
	int curr[N] = { 0 };
	while (true) {
		for (int i = N - 1; curr[i] == 0 && i >= 0; i--) {
			std::cout << '[';
		}
		std::cout << operator[](curr);
		if (curr[N - 1] != _dim[N - 1] - 1) {
			std::cout << ", ";
		}
		int i;
		for (i = N - 1; i >= 0; i--) {
			curr[i]++;
			if (curr[i] == _dim[i]) {
				curr[i] = 0;
				std::cout << ']';
				if (i == 0) {
					std::cout << "\n";
					return;
				}
			} else {
				break;
			}
		}
		if (i < N - 1) {
			std::cout << ",\n";
			if (i < N - 2) {
				std::cout << "\n";
			}
		}
	}
}

template<int N, typename T>
inline void Array<N, T>::save(std::ostream &os) const {
	os << "array-begin" << std::endl;
	os << _initialized << std::endl;
	if (isInitialized()) {
		for (int i = 0; i < N; i++) {
			os << " " << _dim[i];
		}
		os << " " << _size;
		for (int i = 0; i < N; i++) {
			os << " " << offset[i];
		}
		os << std::endl;
		for (auto d : data) {
			os << " " << d;
		}
		os << std::endl;
	}
	os << "array-end" << std::endl;
}

template<int N, typename T>
inline void Array<N, T>::load(std::istream &is) {
	std::string token;
	is >> token;
	if (token != "array-begin") {
		std::cout << "array loading failed" << std::endl;
		exit(0);
	}
	bool init;
	is >> init;
	if (init) {
		_initialized = true;
		int s = 1;
		for (int i = 0; i < N; i++) {
			is >> _dim[i];
			s *= _dim[i];
		}
		is >> _size;
		assert(s == _size);
		for (int i = 0; i < N; i++) {
			is >> offset[i];
		}
		data = std::vector<T>(_size);
		for (int i = 0; i < _size; i++) {
			is >> data[i];
		}
		is >> token;
		if (token != "array-end") {
			std::cout << "unexpected token" << std::endl;
			exit(0);
		}
	} else {
		_initialized = false;
	}
}

/****************************
	DEBUGGING FUNCTIONS
****************************/
template<int N, typename T1, typename T2>
bool check_dim_agree(const Array<N, T1>& t1, const Array<N, T2> &t2) {
	for (int i = 0; i < N; i++) {
		if (t1.dim[i] != t2.dim[i]) {
			return false;
		}
	}
	return true;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v) {
	os << "(";
	for (auto i : v) {
		os << i << ", ";
	}
	os << ")";
	return os;
}