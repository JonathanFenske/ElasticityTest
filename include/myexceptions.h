#ifndef _MYEXCEPTIONS_H_
#define _MYEXCEPTIONS_H_

#include <deal.II/base/exceptions.h>

namespace Elasticity
{
  using namespace dealii;

  // class MyExcDims : public ExceptionBase
  // {
  // public:
  //   MyExcDims(const unsigned int a1,
  //             const unsigned int a2,
  //             const unsigned int a3)
  //     : arg1(a1)
  //     , arg2(a2)
  //     , arg3(a3)
  //   {}
  //   virtual void
  //   print_info(std::ostream &out) const
  //   {
  //     out << "    " outsequence << std::endl;
  //   }

  // private:
  //   unsigned int arg1;
  //   unsigned int arg2;
  //   unsigned int arg3;
  // };

  DeclException3(MyExcDims,
                 unsigned int,
                 unsigned int,
                 unsigned int,
                 "The dimension is " << arg1 << "but must be either " << arg2
                                     << " or " << arg3);
} // namespace Elasticity

#endif // _MYEXCEPTIONS_H_