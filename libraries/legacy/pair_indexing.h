/****************************************************************
  pair_indexing.h                       

  Defines matrix triangle indexing scheme.
                                  
  Mark A. Caprio, University of Notre Dame.

  4/15/11 (mac): Originated.
  4/23/11 (mac): Basic range checking via assertions.
  4/25/15 (mac): Reformatted source file.
  8/26/16 (mac): Minimal patches to serve as deprecated legacy 
    library w/in shell project.

****************************************************************/

#ifndef pair_indexing_h
#define pair_indexing_h

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace legacy {

  ////////////////////////////////////////////////////////////////
  // predefined constants
  ////////////////////////////////////////////////////////////////

  // indexing label

  enum IndexingType {kLowerTriangular, kUpperTriangular, kSquare};

  ////////////////////////////////////////////////////////////////
  // dimension-independent indexing functions
  ////////////////////////////////////////////////////////////////

  // Note: These indexing schemes start from the UL corner and iterate
  // between the diagonal and an "inner edge" (index 0) instead of an
  // "outer edge" (index D-1), leading to a dimension-independent
  // mapping (i,j) -> k.

  // lower triangle
  //   >
  //   ->
  //   -->

  inline int LowerTriangularIndex(int i, int j)
  {
    assert( (0<=j) && (j<=i) );
    return i*(i+1)/2 + j;
  }

  // upper triangle
  //   v||
  //    v|
  //     v

  inline int UpperTriangularIndex(int i, int j)
  {
    assert( (0<=i) && (i<=j) );
    return LowerTriangularIndex(j, i);
  }

  // square
  //   .^^
  //   >.|
  //   ->.  

  inline int SquareIndex(int i, int j)
  {
    assert( (0<=i) && (0<=j) );
    return (j <= i) ? i*i+j : j*(j+2)-i;
  }

  ////////////////////////////////////////////////////////////////
  // dimension-dependent indexing functions -- lexicographical
  ////////////////////////////////////////////////////////////////

  // Note: These indexing schemes traverse (i,j) preserving
  // lexicographical order, at the cost of giving a dimension-dependent
  // indexing.

  // upper triangle
  //   -->
  //    ->
  //     >

  inline int UpperTriangularIndexLex(int D, int i, int j)
  {
    assert( (0<=i) && (i<=j) && (j < D) );
    return i*(2*D-1-i)/2+j;
  }

  // square
  //   -->
  //   -->
  //   -->  

  inline int SquareIndexLex(int D, int i, int j)
  {
    assert( (0<=i) && (0<=j) && (i < D) && (j < D) );
    return i*D+j;
  }

  ////////////////////////////////////////////////////////////////
  // dimension formulas
  ////////////////////////////////////////////////////////////////

  inline int TriangularDimension(int D)
  {
    assert ( (0<=D) );
    return D*(D+1)/2;
  }

  inline int SquareDimension(int D)
  {
    assert ( (0<=D) );
    return D*D;
  }

  ////////////////////////////////////////////////////////////////
  // indexing wrappers
  ////////////////////////////////////////////////////////////////

  inline int LexicographicIndex(IndexingType indexing_type, int D, int i, int j) {
    switch (indexing_type)
      {
      case kLowerTriangular:
	return LowerTriangularIndex(i,j);
	break;
      case kUpperTriangular:
	return UpperTriangularIndexLex(D,i,j);
	break;
      case kSquare:
	return SquareIndexLex(D,i,j);
	break;
      }
  }

  inline int IndexingDimension(IndexingType indexing_type, int D) {
    switch (indexing_type)
      {
      case kLowerTriangular:
      case kUpperTriangular:
	return TriangularDimension(D);
	break;
      case kSquare:
	return SquareDimension(D);
	break;
      }
  }

  ////////////////////////////////////////////////////////////////
  // two-index access to data
  // (a.k.a. a square or triangular array)
  ////////////////////////////////////////////////////////////////

  // Note: coding assumes T is a simple type and therefore passes and returns by value

  template <typename T>
    class PairLookupArray {

  public:

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // default
    // Note: subsequent call to an initializer required for well-defined behavior
  PairLookupArray() : indexing_type_(kSquare), minor_dimension_(0), major_dimension_(0) {};

    // copy -- synthesized

    // initializing
    PairLookupArray (IndexingType indexing_type, int minor_dimension, T value){
      Initialize(indexing_type,minor_dimension,value);
    };

    ////////////////////////////////////////////////////////////////
    // allocation and deallocation
    ////////////////////////////////////////////////////////////////

    // allocate (via resize) and value-initialize
    void Initialize (IndexingType indexing_type, int minor_dimension, T value){
      indexing_type_ = indexing_type;
      minor_dimension_ = minor_dimension;
      major_dimension_ = IndexingDimension(indexing_type, minor_dimension);
      value_vector_.resize(0);
      value_vector_.resize(major_dimension_, value);
    };

    // release allocation
    void Free (){
      indexing_type_ = kSquare;
      minor_dimension_ = 0;
      major_dimension_ = 0;
      value_vector_.clear();
    };

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // entry direct access

    T& operator () (int i, int j) {
      return value_vector_[InternalIndex(i,j)];
    };

    const T& operator () (int i, int j) const {
      return value_vector_[InternalIndex(i,j)];
    };

    // entry get/set access

    T GetValue (int i, int j) const {
      return value_vector_[InternalIndex(i,j)];
    };

    void SetValue (int i, int j, T value) {
      value_vector_[InternalIndex(i,j)]=value;
    };

    ////////////////////////////////////////////////////////////////
    // indexing
    ////////////////////////////////////////////////////////////////

    int MinorDimension () const {return minor_dimension_;};
    int MajorDimension () const {return major_dimension_;};
    int InternalIndex (int i, int j)  const {return LexicographicIndex(indexing_type_,minor_dimension_,i,j);};
    // int DEBUGInternalIndex (int i, int j)  const {
    // 	std::cout << indexing_type_ << " " << minor_dimension_ << " " << i << " " << j << std::endl;
    // 	return LexicographicIndex(indexing_type_,minor_dimension_,i,j);
    // };

  private:

    ////////////////////////////////////////////////////////////////
    // data
    ////////////////////////////////////////////////////////////////

    IndexingType indexing_type_;	
    int minor_dimension_, major_dimension_;
    std::vector< T > value_vector_;

  };



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
