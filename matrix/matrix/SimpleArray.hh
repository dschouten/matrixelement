//
// author: Doug Schouten <doug dot schouten at triumf dot ca>
//

#ifndef SIMPLEARRAY_HH
#define SIMPLEARRAY_HH

template<typename T, unsigned SIZE>
class SimpleArray
{
   public:
      SimpleArray() {}
      explicit SimpleArray(T initial)
      {
         for (unsigned i = 0; i < SIZE; ++i)
         {
            array[i] = initial;
         }
      }

      SimpleArray(const SimpleArray& rhs)
      {
         for (unsigned i = 0; i < SIZE; ++i)
         {
            array[i] = rhs.array[i];
         }
      }

      const SimpleArray& operator=(const SimpleArray& rhs)
      {
         for (unsigned i = 0; i < SIZE; ++i)
         {
            array[i] = rhs.array[i];
         }
         return *this;
      }

      T& operator[](int index)
      {
         return array[index];
      }
      
      const T& operator[](int index) const
      {
         return array[index];
      }

      T array[SIZE];
};

#endif
