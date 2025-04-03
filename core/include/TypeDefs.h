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
typedef ROOT::Math::SMatrix<double, 1, 5> SMatrix15;
typedef ROOT::Math::SMatrix<double, 2, 5> SMatrix25;
typedef ROOT::Math::SMatrix<double, 2, 7> SMatrix27;
typedef ROOT::Math::SMatrix<double, 2, 12> SMatrix212;
typedef ROOT::Math::SMatrix<double, 3, 2> SMatrix32;
typedef ROOT::Math::SMatrix<double, 3, 3> SMatrix33;
typedef ROOT::Math::SMatrix<double, 3, 6> SMatrix36;
typedef ROOT::Math::SMatrix<double, 5, 1> SMatrix51;

typedef ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double, 2> >    SMatrixSym2;
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> >    SMatrixSym3;
typedef ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5> >    SMatrixSym5;
typedef ROOT::Math::SMatrix<double, 6, 6, ROOT::Math::MatRepSym<double, 6> >    SMatrixSym6;

typedef ROOT::Math::SVector<double, 1> SVector1;
typedef ROOT::Math::SVector<double, 2> SVector2;
typedef ROOT::Math::SVector<double, 5> SVector5;
typedef ROOT::Math::SVector<double, 6> SVector6;