#pragma once

#include<iostream>
#include<assert.h>
#include<stdexcept>
#include<vector>
#include<unordered_map>


namespace me{
	template<typename T>
	inline void swap(T& a, T& b) {
		T hold = a;
		a = b;
		b = hold;
		return;
	}

	template<typename T>
	inline bool is_zero(T& a) {
		if (a > -10E-6 && a < 10E-6) {
			return true;
		}
		else {
			return false;
		}
	}

	inline void handle_memory_alloc(void* point_to_allocated_memory) {
		assert(!(point_to_allocated_memory == nullptr));
	}

	template<typename T>
	inline constexpr bool is_decimal(const double& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_decimal(const float& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_decimal(const long double& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_decimal(const T& v)
	{
		return false;
	}

	template<typename T>
	inline constexpr bool is_integer(const int& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_integer(const short& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_integer(const long& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_integer(const char& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_integer(const long long& v)
	{
		return true;
	}
	template<typename T>
	inline constexpr bool is_integer(const T& v)
	{
		return false;
	}
}
