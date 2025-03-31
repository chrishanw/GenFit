/*
 * Genfit typedefs
 * TypeDefs.h
 *
 *  Created on: Mar 30, 2025
 *      Author: C. Wessel
 */

#include <Math/SMatrix.h>
#include <Math/SVector.h>
typedef ROOT::Math::SMatrix<double, 2> SMatrix22;
typedef ROOT::Math::SMatrix<double, 3> SMatrix33;
typedef ROOT::Math::SMatrix<double, 5> SMatrix55;
typedef ROOT::Math::SMatrix<double, 2, 3> SMatrix23;
typedef ROOT::Math::SMatrix<double, 2, 5> SMatrix25;
typedef ROOT::Math::SMatrix<double, 2, 7> SMatrix27;
typedef ROOT::Math::SMatrix<double, 3, 2> SMatrix32;
typedef ROOT::Math::SMatrix<double, 3, 6> SMatrix36;

typedef ROOT::Math::SVector<double, 2> SVector2;
typedef ROOT::Math::SVector<double, 5> SVector5;

typedef ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double, 2> >    SMatrixSym2;
typedef ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5> >    SMatrixSym5;