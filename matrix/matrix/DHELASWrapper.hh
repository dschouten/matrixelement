#ifndef DHELASWRAPPER_HH
#define DHELASWRAPPER_HH

#include <complex>
#include <vector>
#include <iostream>
#include <TLorentzVector.h>

#include "matrix/SimpleArray.hh"
#include "matrix/Utility.hh"

typedef std::complex<double> doublecomplex;

extern "C"
{
   int ixxxxx_(double*, double*, int*, int*, doublecomplex*);
   int oxxxxx_(double*, double*, int*, int*, doublecomplex*);
   int vxxxxx_(double*, double*, int*, int*, doublecomplex*);

   int jioxxx_(doublecomplex*, doublecomplex*, doublecomplex*, double*,
               double*, doublecomplex*);
   int fvoxxx_(doublecomplex*, doublecomplex*, doublecomplex*, double*,
               double*, doublecomplex*);
   int fvixxx_(doublecomplex*, doublecomplex*, doublecomplex*, double*,
               double*, doublecomplex*);
   int jvvxxx_(doublecomplex*, doublecomplex*, doublecomplex*, double*,
               double*, doublecomplex*);
   int hvvxxx_(doublecomplex*, doublecomplex*, doublecomplex*, double*,
               double*, doublecomplex*);
   int hioxxx_(doublecomplex*, doublecomplex*, doublecomplex*, double*,
               double*, doublecomplex*);
   int jvsxxx_(doublecomplex*, doublecomplex*, doublecomplex*, double*,
               double*, doublecomplex*);

   int iovxxx_(doublecomplex*, doublecomplex*, doublecomplex*, doublecomplex*,
               doublecomplex*);
   int iosxxx_(doublecomplex*, doublecomplex*, doublecomplex*, doublecomplex*,
               doublecomplex*);
   int vvvxxx_(doublecomplex*, doublecomplex*, doublecomplex*, doublecomplex*,
               doublecomplex*);
}

namespace DHELAS
{
   struct HelArray
   {
         HelArray() {p1 = p2 = p3 = p4 = p5 = p6 = 0;}
         doublecomplex p1, p2, p3, p4, p5, p6;
   };

   template<unsigned HELICITY>
   inline SimpleArray<HelArray, HELICITY>
   inputFunction(TLorentzVector input, double mass, int direction, 
                 int(*function)(double*, double*, int*, int*, doublecomplex*),
                 bool rh = false);

   template<>
   inline SimpleArray<HelArray, 1>
   inputFunction<1>(TLorentzVector input, double mass, int direction,
                    int(*function)(double*, double*, int*, int*,
                                   doublecomplex*),
                    bool rh);
   template<>
   inline SimpleArray<HelArray, 2>
   inputFunction<2>(TLorentzVector input, double mass, int direction,
                    int(*function)(double*, double*, int*, int*,
                                   doublecomplex*),
                    bool rh);
   template<>
   inline SimpleArray<HelArray, 3>
   inputFunction<3>(TLorentzVector input, double mass, int direction,
                    int(*function)(double*, double*, int*, int*,
                                   doublecomplex*),
                    bool rh);

   template<unsigned SIZE1, unsigned SIZE2>
   inline SimpleArray<HelArray, SIZE1 * SIZE2>
   middleFunction(SimpleArray<HelArray, SIZE1>& input1,
                  SimpleArray<HelArray, SIZE2>& input2,
                  doublecomplex factor[2], double mass, double width,
                  int(*function)(doublecomplex*, doublecomplex*,
                                 doublecomplex*, double*, double*,
                                 doublecomplex*));

   template<unsigned SIZE1, unsigned SIZE2, unsigned SIZE3>
   inline SimpleArray<doublecomplex, SIZE1 * SIZE2 * SIZE3>
   outputFunction(SimpleArray<HelArray, SIZE1>& input1,
                  SimpleArray<HelArray, SIZE2>& input2,
                  SimpleArray<HelArray, SIZE3>& input3,
                  doublecomplex factor[2],
                  int(*function)(doublecomplex*, doublecomplex*,
                                 doublecomplex*, doublecomplex*,
                                 doublecomplex*));

   template<unsigned HELICITY>
   SimpleArray<HelArray, HELICITY>
   ixxxxx(TLorentzVector input, double mass, int direction, bool rh = false)
   {
      return inputFunction<HELICITY>(input, mass, direction, ixxxxx_, rh);
   }
   template<unsigned HELICITY>
   SimpleArray<HelArray, HELICITY>
   oxxxxx(TLorentzVector input, double mass, int direction, bool rh = false)
   {
      return inputFunction<HELICITY>(input, mass, direction, oxxxxx_, rh);
   }
   template<unsigned HELICITY>
   SimpleArray<HelArray, HELICITY>
   vxxxxx(TLorentzVector input, double mass, int direction)
   {
      return inputFunction<HELICITY>(input, mass, direction, vxxxxx_);
   }

   template<unsigned SIZE1, unsigned SIZE2>
   SimpleArray<HelArray, SIZE1 * SIZE2>
   jioxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          doublecomplex factor[2], double mass, double width)
   {
      return middleFunction(input1, input2, factor, mass, width, jioxxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2>
   SimpleArray<HelArray, SIZE1 * SIZE2>
   fvoxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          doublecomplex factor[2], double mass, double width)
   {
      return middleFunction(input1, input2, factor, mass, width, fvoxxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2>
   SimpleArray<HelArray, SIZE1 * SIZE2>
   fvixxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          doublecomplex factor[2], double mass, double width)
   {
      return middleFunction(input1, input2, factor, mass, width, fvixxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2>
   SimpleArray<HelArray, SIZE1 * SIZE2>
   jvvxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          doublecomplex factor[2], double mass, double width)
   {
      return middleFunction(input1, input2, factor, mass, width, jvvxxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2>
   SimpleArray<HelArray, SIZE1 * SIZE2>
   hvvxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          doublecomplex factor[2], double mass, double width)
   {
      return middleFunction(input1, input2, factor, mass, width, hvvxxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2>
   SimpleArray<HelArray, SIZE1 * SIZE2>
   hioxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          doublecomplex factor[2], double mass, double width)
   {
      return middleFunction(input1, input2, factor, mass, width, hioxxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2>
   SimpleArray<HelArray, SIZE1 * SIZE2>
   jvsxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          doublecomplex factor[2], double mass, double width)
   {
      return middleFunction(input1, input2, factor, mass, width, jvsxxx_);
   }

   template<unsigned SIZE1, unsigned SIZE2, unsigned SIZE3>
   SimpleArray<doublecomplex, SIZE1 * SIZE2 * SIZE3>
   iovxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          SimpleArray<HelArray, SIZE3>& input3,
          doublecomplex factor[2])
   {
      return outputFunction(input1, input2, input3, factor, iovxxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2, unsigned SIZE3>
   SimpleArray<doublecomplex, SIZE1 * SIZE2 * SIZE3>
   iosxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          SimpleArray<HelArray, SIZE3>& input3,
          doublecomplex factor[2])
   {
      return outputFunction(input1, input2, input3, factor, iosxxx_);
   }
   template<unsigned SIZE1, unsigned SIZE2, unsigned SIZE3>
   SimpleArray<doublecomplex, SIZE1 * SIZE2 * SIZE3>
   vvvxxx(SimpleArray<HelArray, SIZE1>& input1,
          SimpleArray<HelArray, SIZE2>& input2,
          SimpleArray<HelArray, SIZE3>& input3,
          doublecomplex factor[2])
   {
      return outputFunction(input1, input2, input3, factor, vvvxxx_);
   }
}

std::ostream& operator<<(std::ostream& stream, const DHELAS::HelArray& rhs);

#include "DHELASWrapper.icc"

#endif
