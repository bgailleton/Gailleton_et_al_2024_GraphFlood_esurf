/*
Universal vector-like data holder, used to make the data type generic.
It simplifies the original approach of always calling format_input/format_output
at each functions These will only be called here and all the data holders in the
code will be of that type.
*/

#ifndef VECLIKE_HPP
#define VECLIKE_HPP

template<class T>
class veclike
{
public:
	T* ptr = nullptr;
	int isize = 0;
	size_t usize = 0;

	veclike(){};

	template<class U>
	veclike(std::vector<U>& data)
	{
	}

#ifdef DAGGER_FT_PYTHON
	veclike(py::array_t<T, 1>& arr)
	{
		auto boeuf_heure = arr.request();
		this->ptr = (T*)boeuf_heure.ptr;
		this->isize = arr.size();
		this->usize = arr.size();
	};

#endif

	T& operator[](int i) { return this->ptr[i]; }

	void set(int i, T val) { (*this)[i] = val; }
	T get(int i) { return (*this)[i]; }

	std::vector<T> to_vec()
	{
		std::vector<T> out(this->isize);
		for (int i = 0; i < this->isize; ++i)
			out[i] = (*this)[i];
	}

	// size_t size(){return this->usize;}
	const size_t size() { return this->usize; }
};

#endif
