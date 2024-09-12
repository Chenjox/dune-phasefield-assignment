#ifndef BIW4_67_PHASEFIELD_HH
#define BIW4_67_PHASEFIELD_HH


#include <dune/common/fvector.hh> 
#include <dune/common/fmatrix.hh>
// add your classes here


namespace Dune {
    namespace Biw467 {
        

        template<int dim, class field_type=double>
        class material {

            public:

            double _griffithReleaseRate,_shearModulus,_firstLameConstant,_regularisationParameter;

            material(double griffithReleaseRate,double shearModulus,double firstLameConstant,double regularisationParameter) :
                _griffithReleaseRate(griffithReleaseRate),
                _shearModulus(shearModulus),
                _firstLameConstant(firstLameConstant),
                _regularisationParameter(regularisationParameter) {};

            double phasefieldEnergyDensity(double phasefieldNumber, const Dune::FieldVector<field_type, dim> gradient) const {
                double normgrad = gradient.two_norm2();
                return _griffithReleaseRate * 1.0/(2.0 * _regularisationParameter)*(phasefieldNumber*phasefieldNumber + _regularisationParameter*_regularisationParameter*normgrad);
            }

            double degradationFunction(double phasefieldNumber) const {
                return (1.0 - phasefieldNumber)*(1.0 - phasefieldNumber);
            }

            double degradationFunctionDerivative(double phasefieldNumber) const {
                return -2.0*(1.0 - phasefieldNumber);
            }

            double strainEnergyDensity(const Dune::FieldMatrix<field_type, dim, dim>& strains) const {
                
                double trace = 0.0;
                for (int i = 0; i <dim; i++) {
                    trace += strains[i][i];
                }
                if (dim == 2) { // Plane Strain correction
                    trace += 1.0;
                }
                double traceE2 = strains.frobenius_norm2();

                return _shearModulus*traceE2 + 0.5*_firstLameConstant*trace*trace;
            }

            double stresses(const Dune::FieldMatrix<field_type, dim, dim>& strains, Dune::FieldMatrix<field_type, dim, dim>& stresses) const {
                stresses = 0.0;
                double trace = 0.0;
                for (int i = 0; i <dim; i++) {
                    trace += strains[i][i];
                }
                if (dim == 2) { // Plane Strain correction
                    trace += 1.0;
                }
                for (int i; i < dim; i++) {
                    for (int j; j < dim; j++) {
                        stresses[i][j] += _firstLameConstant * strains[i][j];
                    }
                    stresses[i][i] += 2.0 * _shearModulus * trace;
                }
            }

        };
    }
}


#endif // BIW4_67_PHASEFIELD_HH
