/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_AbsHMatrix_h
#define genfit_AbsHMatrix_h

#include <SMatrixTypeDefs.h>
#include <Math/MatrixFunctions.h>

#include <TMatrixDSym.h>
#include <TVectorD.h>


namespace genfit {

/**
 * @brief HMatrix for projecting from AbsTrackRep parameters to measured parameters in a DetPlane.
 *
 */
template<unsigned int nRows>
class AbsHMatrix {

  using SMatrixX5 = ROOT::Math::SMatrix<double, nRows, 5>;
  using SMatrix5X = ROOT::Math::SMatrix<double, 5, nRows>;
  using SVectorX = ROOT::Math::SVector<double, nRows>;
  using SMatrixSymXD = ROOT::Math::SMatrix<double, nRows, nRows, ROOT::Math::MatRepSym<double, nRows> >;

 public:

  AbsHMatrix() {;}

  virtual ~AbsHMatrix() {;}

  //! Get the actual matrix representation
  virtual const SMatrixX5& getMatrix() const = 0;

  //! H*v
  virtual SVectorX Hv(const SVector5& v) const {return getMatrix()*v;}

  //! M*H^t
  virtual SMatrix5X MHt(const SMatrixSym5& M) const {
    // return TMatrixD(M, TMatrixD::kMultTranspose, getMatrix());
    return M * ROOT::Math::Transpose(getMatrix());
  }

  template<unsigned int mRows>
  ROOT::Math::SMatrix<double, mRows, nRows> MHt(const ROOT::Math::SMatrix<double, mRows, 5>& M) const
  {
    // return TMatrixD(M, TMatrixD::kMultTranspose, getMatrix());
    return M * ROOT::Math::Transpose(getMatrix());
  }

  //! similarity: H*M*H^t
  virtual SMatrixSymXD HMHt(SMatrixSym5& M) const {return ROOT::Math::Similarity(M, getMatrix());}

  virtual AbsHMatrix<nRows>* clone() const = 0;

  bool operator==(const AbsHMatrix& other) const {return this->isEqual(other);}
  bool operator!=(const AbsHMatrix& other) const {return !(this->isEqual(other));}
  virtual bool isEqual(const AbsHMatrix& other) const = 0;

  virtual void Print(const Option_t* = "") const {;}

 protected:
  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsHMatrix& operator=(const AbsHMatrix&);

 public:
  ClassDef(AbsHMatrix,2)
  // Version history:
  // ver 2: no longer derives from TObject

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsHMatrix_h
