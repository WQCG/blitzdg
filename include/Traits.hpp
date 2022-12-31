// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Traits.hpp
 * @brief Defines a basic set of traits for making compile-time assertions
 * about numeric types and iterators.
 */
#pragma once
#include <iterator>
#include <type_traits>

namespace blitzdg {
	/////////////////////////////////////////////////////////////////////
	// GENERAL
	/////////////////////////////////////////////////////////////////////
    /**
     * Returns true if T and U are the same type.
     */
    template <typename T, typename U>
    constexpr bool isSame() {
        return std::is_same<T, U>::value;
    }

    /**
     * Returns true if T is an integral numeric type.
     */
    template <typename T>
    constexpr bool isIntegral() {
        return std::is_integral<T>::value;
    }
	
	/**
     * Returns true if T is a floating-point numeric type.
     */
	template <typename T>
	constexpr bool isFloatingPoint() {
		return std::is_floating_point<T>::value;
	}
	
	/**
	 * Returns true if T is either an integral or floating-point
	 * numeric type.
	 */
	template <typename T>
	constexpr bool isArithmetic() {
		return isIntegral<T>() || isFloatingPoint<T>();
	}

	/////////////////////////////////////////////////////////////////////
	// ITERATORS
	/////////////////////////////////////////////////////////////////////
	/**
	 * Defines an alias, type, for the value_type of an iterator.
	 */
	template <typename Iterator>
	struct IteratorValueType {
		using type = typename std::iterator_traits<Iterator>::value_type;
	};

	/**
	 * Specialization of IteratorValueType for std::back_insert_iterator.
	 */
	template <typename Container>
	struct IteratorValueType<std::back_insert_iterator<Container>> {
		using type = typename Container::value_type;
	};

	/**
	 * Returns true if Iterator1 and Iterator2 have the same value_type.
	 */
	template <typename Iterator1, typename Iterator2>
	constexpr bool iteratorValueTypesSame() {
		return std::is_same<
			typename IteratorValueType<Iterator1>::type,
			typename IteratorValueType<Iterator2>::type>::value;
	}

	/**
	 * Returns true if the value_type of Iterator is the as the type T.
	 */
	template <typename Iterator, typename T>
	constexpr bool iteratorValueTypeSameAs() {
		return std::is_same<T,
			typename IteratorValueType<Iterator>::type>::value;
	}

	/**
	 * Returns true if Iterator is a random-access iterator.
	 */
	template <typename Iterator>
	constexpr bool randomAccessIterator() {
		return std::is_same<std::random_access_iterator_tag,
			typename std::iterator_traits<Iterator>::iterator_category>::value;
	}

	/**
	 * Returns true if Iterator is a bidirectional iterator.
	 */
	template <typename Iterator>
	constexpr bool bidirectionalIterator() {
		return std::is_same<std::bidirectional_iterator_tag,
			typename std::iterator_traits<Iterator>::iterator_category>::value ||
			randomAccessIterator<Iterator>();
	}

	/**
	 * Returns true if Iterator is a forward iterator.
	 */
	template <typename Iterator>
	constexpr bool forwardIterator() {
		return std::is_same<std::forward_iterator_tag,
			typename std::iterator_traits<Iterator>::iterator_category>::value ||
			bidirectionalIterator<Iterator>();
	}

	/**
	 * Returns true if Iterator is an input iterator.
	 */
	template <typename Iterator>
	constexpr bool inputIterator() {
		return std::is_same<std::input_iterator_tag,
			typename std::iterator_traits<Iterator>::iterator_category>::value ||
			forwardIterator<Iterator>();
	}

	/**
	 * Returns true if Iterator is an output iterator.
	 */
	template <typename Iterator>
	constexpr bool outputIterator() {
		return std::is_same<std::output_iterator_tag,
			typename std::iterator_traits<Iterator>::iterator_category>::value ||
			forwardIterator<Iterator>();
	}
} // namespace blitzdg